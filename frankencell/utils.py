
from itertools import count
import h5py as h5
import pandas as pd
import numpy as np
import anndata
from scipy import sparse
import tempfile
import os
from scipy.io import mmwrite
import subprocess
from dynclipy.read import ro as dynro
from rpy2 import rinterface

@dynro.conversion.rpy2py.register(rinterface.NULLType)
def convert_NULL(obj):
    return None


def decode_h5(v):
    if type(v[0]) == np.bytes_:
        return np.array([x.decode() for x in v])
    else:
        return v

def read_info(dir):

    columns = list(dir.keys())
    info = pd.DataFrame(
        {col : decode_h5(dir[col][...]) for col in columns}
    )

    return info

def read_cell_info(filename):

    with h5.File(filename, 'a') as f:
        
        dir = f['data']['cell_info']['data']
        return read_info(dir)

def read_feature_info(filename):

    with h5.File(filename, 'a') as f:
        
        dir = f['data']['feature_info']['data']
        return read_info(dir)


def dataset_to_anndata(dataset):
    
    counts = sparse.csr_matrix(dataset['counts'].values).astype(np.float32)

    return anndata.AnnData(
        obs = dataset['cell_info'].set_index('cell_id'),
        var = dataset['feature_info'].set_index('feature_id'),
        obsm = {'dimred' : dataset['dimred']} if 'dimred' in dataset.keys() else None,
        layers = {
            'counts' : counts,
        },
        X = dataset['expression'].values if not dataset['expression'] is None else counts
    )

def read_dynframe(dynframe_path):

    dataset = dynro.r('dynutils::read_h5("{}")'.format(dynframe_path))
    return dataset_to_anndata(dataset)


def select_features(adata, feature_type):

    if 'feature_type'  in adata.var.columns:
        if not ((adata.var_vector('feature_type') == feature_type).sum() > 0):
            raise ValueError('Feature type {} not in dataset.'.format(feature_type))

        adata = adata[:, adata.var_vector('feature_type') == feature_type]

    else:
        raise IndexError('Adata not labeled with feature types.')

    return adata

def add_expression_to_dynframe(dynframe, output_file, feature_df, expression = None, 
        counts = None,
    ):

    assert(not expression is None or not counts is None)
    with tempfile.TemporaryDirectory() as folder:
        
        if not counts is None:
            assert(sparse.isspmatrix(counts))
            counts_file = os.path.join(folder, 'counts')
            mmwrite(counts_file, counts)
            counts_file+='.mtx'
        else:
            counts_file = None
        
        if not expression is None:
            if isinstance(expression, np.ndarray) or not sparse.issparse(expression):
                expression = sparse.csr_matrix(expression)

            expr_file = os.path.join(folder, 'expr.csv')
            mmwrite(expr_file, expression)
            expr_file+='.mtx'
        else:
            expr_file = None
        
        featureinfo_file = os.path.join(folder, 'feature_info.csv')
        feature_df.to_csv(featureinfo_file, index = None)

        rstring = '''
dataset <- dynutils::read_h5("{dynframe}")
feature_info = tibble::tibble(
    read.csv("{featureinfo_file}")
)

if ({use_expr_file}) {{
    expr <- Matrix::readMM("{expr_file}")
    rownames(expr) <- dataset$cell_ids
    colnames(expr) <- feature_info$feature_id
}} else {{
    expr <- NULL
}}


if ({use_counts_file}) {{
    counts <- Matrix::readMM("{counts_file}")
    rownames(counts) <- dataset$cell_ids
    colnames(counts) <- feature_info$feature_id
}} else {{
    counts <- NULL
}}

dataset <- dynwrap::add_expression(dataset, counts = counts, expression = expr,
                                   feature_info = feature_info)

dynutils::write_h5(dataset, "{output_file}")

        '''.format(
            dynframe = dynframe, counts_file = counts_file,
            expr_file = expr_file, featureinfo_file = featureinfo_file,
            output_file = output_file, 
            use_expr_file = str(not expr_file is None).upper(),
            use_counts_file = str(not counts_file is None).upper(),
        )

        dynro.r(rstring)


def write_cell_info(filename, info_dict):

    info_dict = {
        k : v.astype('S') if np.issubdtype(v.dtype, str) else v
        for k, v in info_dict.items()
    }

    with h5.File(filename, 'a') as f:

        colnames = np.array(['cell_id', *info_dict.keys()]).astype('S')
        
        dir = f['data']['cell_info']
        del dir['colnames']
        dir.create_dataset("colnames", data = colnames)

        for key, val in info_dict.items():
            if key in dir['data']:
                del dir['data'][key]

            dir['data'].create_dataset(key, data = val)


def add_dimred_prior(filename, dimred):

    with tempfile.TemporaryDirectory() as folder:
                
        dimred_path = os.path.join(folder, 'dimred.csv')
        np.savetxt(dimred_path, dimred)
        
        rstring='''
dataset <- dynutils::read_h5("{dynpath}")
dimred <- as.matrix(
    read.table("{dimred_path}", sep = ' ', header = FALSE)
)
rownames(dimred) <- dataset$cell_ids

dataset <- dynwrap::add_prior_information(dataset, dimred = dimred)
dynutils::write_h5(dataset, "{dynpath}")
        '''.format(dynpath = filename, dimred_path = dimred_path)
        
        dynro.r(rstring)
