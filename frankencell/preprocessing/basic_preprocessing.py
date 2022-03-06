import scanpy as sc

def basic_rna_preprocessing(adata, min_cells, min_dispersion):

    sc.pp.filter_genes(adata, min_cells = min_cells)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_disp=min_dispersion)

    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata)

    return adata


def basic_atac_preprocessing(adata, min_cells):
    sc.pp.filter_genes(adata, min_cells = min_cells)

    return adata