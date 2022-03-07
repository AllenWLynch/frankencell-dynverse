
from ..utils import read_dynframe, select_features, add_expression_to_dynframe
from .basic_preprocessing import basic_atac_preprocessing
import logging
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.decomposition import TruncatedSVD

def main(
    dynframe_path,
    output_path,
    feature_type = 'ATAC',
    min_cells = 25,
):

    adata = read_dynframe(dynframe_path)

    try:
        adata = select_features(adata, feature_type)
    except IndexError:
        logging.warning('Cannot find feature type {}, assuming all features are {}.'.format(feature_type, feature_type))

    adata = basic_atac_preprocessing(adata, min_cells)

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
    parser.add_argument('--feature-type', '-f', type = str, default = 'ATAC')
    parser.add_argument('--min-cells', '-m', type = int, default = 25)