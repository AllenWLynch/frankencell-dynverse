
from ..utils import read_dynframe, select_features, add_expression_to_dynframe
from .basic_preprocessing import basic_atac_preprocessing, basic_rna_preprocessing
import anndata
from scipy import sparse
import numpy as np

def main(
    dynframe_path,
    output_path,
    min_cells = 25,
    min_dispersion = 0.7
):

    adata = read_dynframe(dynframe_path)

    atac_data = select_features(adata, 'ATAC')
    atac_data = basic_atac_preprocessing(atac_data, min_cells)

    rna_data = select_features(adata, 'RNA')
    rna_data = basic_rna_preprocessing(rna_data, min_cells, min_dispersion=min_dispersion)

    divcol = sparse.csr_matrix(np.zeros((len(rna_data), 1)))
    divider = anndata.AnnData(
        X = divcol, layers = {'counts' : divcol},
        obs = rna_data.obs    
    )

    adata = anndata.concat([rna_data, divider, atac_data], axis = 1)

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        adata.var.reset_index(), 
        counts = adata.layers['counts'], 
        expression = adata.X.copy()
    )

    return adata

def add_arguments(parser):

    parser.add_argument('--dynframe', '-i', type = str, required = True)
    parser.add_argument('--outfile', '-o', type = str, required = True)
    parser.add_argument('--feature-type', '-f', type = str, default = 'ATAC')
    parser.add_argument('--min-cells', '-m', type = int, default = 25)