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