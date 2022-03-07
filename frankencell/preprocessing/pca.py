import argparse
import scanpy as sc
from ..utils import read_dynframe, select_features, add_expression_to_dynframe
from .basic_preprocessing import basic_rna_preprocessing
import logging

def main(
    dynframe_path,
    output_path,
    feature_type = 'RNA',
    min_cells = 25,
    min_dispersion = 0.7
):

    adata = read_dynframe(dynframe_path)

    try:
        adata = select_features(adata, feature_type)
    except IndexError:
        logging.warning('Cannot find feature type {}, assuming all features are {}.'.format(feature_type, feature_type))

    adata = basic_rna_preprocessing(adata, min_cells, min_dispersion)

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        adata.var.reset_index(), 
        counts = adata.layers['counts'],
        expression = adata.X.copy()
    )

def add_arguments(parser):
    
    parser.add_argument('--dynframe', '-i', type = str, required = True)
    parser.add_argument('--outfile', '-o', type = str, required = True)
    parser.add_argument('--feature-type', '-f', type = str, default = 'RNA')
    parser.add_argument('--min-cells', '-m', type = int, default = 25)
    parser.add_argument('--min-dispersion', '-d', type = float, default = 0.7)
