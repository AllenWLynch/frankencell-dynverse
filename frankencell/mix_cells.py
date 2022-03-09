from ast import parse
from typing_extensions import Required
from scipy.stats import multivariate_hypergeom, lognorm
from scipy import sparse
import numpy as np
import anndata
from functools import partial
from joblib import Parallel, delayed
from tqdm import tqdm
from .utils import read_cell_info, add_expression_to_dynframe
import pandas as pd
import shutil

from dynclipy.dataset import add_adder, Dataset
add_adder(Dataset, "add_expression")

def prepare_adata(adata, pure_states, 
                  cell_state_col = 'leiden',
                 counts_layer = 'counts'):

    cluster_counts = [
        adata[adata.obs[cell_state_col] == cluster].layers[counts_layer].tocsr()
        for cluster in pure_states
    ]

    cluster_rd = [
        np.array(counts.sum(-1)).reshape(-1)
        for counts in cluster_counts
    ]
    
    return cluster_counts, cluster_rd


def mix_cell_reads(read_depth, mixing_weights, cells):
        
    def geom_sample_sparse_array(arr, n_samples):

        subsampled_reads = multivariate_hypergeom(arr.data, n_samples).rvs()
        return sparse.csr_matrix((subsampled_reads.reshape(-1), arr.indices, arr.indptr), shape = arr.shape)
        
    mixed_cell = sparse.csr_matrix(sparse.vstack([
        geom_sample_sparse_array(feature_counts, int(read_depth * m))
        for feature_counts, m in zip(map(lambda x : x.tocsr().astype(int), cells), mixing_weights)
    ]).sum(0))
    
    return mixed_cell


def get_valid_cells_dataset(
    cluster_rd, mixing_weights, read_depth
):
    return [
        cell_read_depths >= int(read_depth * mix_proportion)
        for cell_read_depths, mix_proportion in zip(cluster_rd, mixing_weights)
    ] # K, M


def get_mixable_cells(cluster_rds, read_depths, mixing_weights, rd_means, rd_stds):

    #D, K, M
    num_clusters = mixing_weights.shape[0]
    num_modes = len(cluster_rds)

    found_valid_cells = False

    while not found_valid_cells:

        valid_cells = list(map(
                lambda x : get_valid_cells_dataset(x[0], mixing_weights, x[1]), 
                zip(cluster_rds, read_depths)
            )) # D, K, M

        valid_cells_per_cluster = [] # K, M

        for cluster in range(num_clusters):
            
            valid_cells_per_cluster.append(
                np.hstack([
                    valid_cells[mode][cluster][:, np.newaxis]
                    for mode in range(num_modes)
                ]).all(-1)
            )

        mode_has_valid_cells = []
        for mode in range(num_modes):
            mode_has_valid_cells.append(
                np.all([valid_cells[mode][cluster].any()
                    for cluster in range(num_clusters)])
            )

        found_valid_cells = np.all([valid_cells.sum() > 0 for valid_cells in valid_cells_per_cluster])

        if not found_valid_cells:
            read_depths = [
                rd if has_valid_cells else sample_bounded_read_depth(rd_mean, rd_std, 2.5, 10)[0]
                for rd, has_valid_cells, rd_mean, rd_std in zip(
                    read_depths, mode_has_valid_cells, rd_means, rd_stds
                )
            ]

    mix_cells = [
        np.random.choice(np.argwhere(valid_cells)[:,0])
        for valid_cells in valid_cells_per_cluster
    ] #K

    return mix_cells, read_depths


def mix_reads(*,cluster_counts, read_depths, cells, mixing_weights): #D, K, M and K

    dataset_mixology_reads = []
    for dataset_clusters, read_depth in zip(cluster_counts, read_depths): #D

        reads = [
            dataset_clusters[clusternum][cellnum]
            for clusternum, cellnum in enumerate(cells) #
        ] # K, G

        mixed_reads = mix_cell_reads(read_depth, mixing_weights, reads) # G

        dataset_mixology_reads.append(
            mixed_reads
        )

    return dataset_mixology_reads #D


def sample_proportions(
    read_depths, mixing_weights,*,
    cluster_counts, 
    cluster_rds,
    rd_means, rd_stds):

    cells, read_depths = get_mixable_cells(cluster_rds, read_depths, mixing_weights, rd_means, rd_stds)

    return mix_reads(
        cluster_counts = cluster_counts, 
        read_depths = read_depths, 
        mixing_weights = mixing_weights,
        cells = cells
    )


def get_bounded_counts(exp_mean, std, max_std, lower = True):
    return np.exp(
        np.log(exp_mean) + (-1 if lower else 1) * std * max_std
    )


def sample_bounded_read_depth(exp_mean, std, max_std, n_cells):

    min_counts, max_counts = get_bounded_counts(exp_mean, std, max_std, lower=True),\
            get_bounded_counts(exp_mean, std, max_std, lower=False)

    found_depths = False
    mult = 2
    while not found_depths:
        samples = lognorm(std, 0, exp_mean).rvs(n_cells * mult)
        allowed_samples = (samples < max_counts) & (samples > min_counts)

        if allowed_samples.sum() >= n_cells:
            read_depths = samples[np.random.choice(np.argwhere(allowed_samples)[:,0], size = n_cells)]
            return read_depths
        else:
            mult +=1

        if mult > 10:
            raise ValueError("Cannot sample enough read depths to meet your 'max_std' requirement. Loosen this requirement or decrease the std of the read depth distribution.")
        
def read_datasets(paths):
    datasets = [anndata.read_h5ad(path) for path in paths]

    return datasets

def check_datasets(datasets, cell_state_col):

    assert( #same length
        np.all(np.array([len(d) for d in datasets]) == len(datasets[0]))
    )

    assert( #same barcodes
        all(
            [np.all(d.obs_names.values == datasets[0].obs_names.values)
            for d in datasets]
        )
    )

    assert( #same cluster column
        all(
            [np.all(d.obs_vector(cell_state_col).astype(str) == datasets[0].obs_vector(cell_state_col).astype(str))
            for d in datasets]
        )
    )


def mix_frankencells(*,
    scaffold,
    datasets, 
    rd_means, 
    rd_stds, 
    pure_states, 
    feature_types,
    max_std = 2.5, 
    cell_state_col = 'leiden',
    counts_layer = 'counts',
    n_jobs = 1,
    seed = None,
    output_path = None,
):

    if output_path is None:
        output_path = scaffold

    shutil.copyfile(scaffold, output_path)

    cell_info = read_cell_info(output_path)
    mixing_weights = cell_info.iloc[:, cell_info.columns.str.startswith('mix_weight_')].values

    assert(np.all([
            len(arg) == len(datasets) 
            for arg in [rd_means, rd_stds, feature_types]]))

    assert(len(pure_states) >= mixing_weights.shape[1])
    datasets = read_datasets(datasets)
    check_datasets(datasets, cell_state_col)

    np.random.seed(seed)

    required_read_depths = list(zip(*[
        sample_bounded_read_depth(mean, std, max_std, len(mixing_weights))
        for mean, std in zip(rd_means, rd_stds)
    ]))

    cluster_counts, cluster_rds = list(zip(*map(lambda d : prepare_adata(d, pure_states, 
                  cell_state_col = cell_state_col,
                 counts_layer = counts_layer), datasets)
            ))
    
    franken_function = partial(sample_proportions,
            cluster_counts = cluster_counts, cluster_rds = cluster_rds,
            rd_means = rd_means, rd_stds = rd_stds)
    
    if n_jobs > 1:
        frankencells = Parallel(n_jobs=n_jobs, verbose = 0, pre_dispatch='2*n_jobs')(
            delayed(franken_function)\
                (mixing_weights = weights, read_depths = rds)
                for weights, rds in tqdm(
                        zip(mixing_weights, required_read_depths), 
                        total = len(mixing_weights), desc = 'Stitching cells'
                )
            )
    else:
        frankencells = [
            franken_function(mixing_weights = weights, read_depths = rds)
            for weights, rds in tqdm(
                zip(mixing_weights, required_read_depths), 
                total = len(mixing_weights),  desc = 'Stitching cells')
        ]
        
    counts_by_mode = list(map(sparse.vstack, list(zip(*frankencells))))

    all_counts = sparse.hstack(counts_by_mode)

    feature_type = [
        x for feature, matrix in zip(feature_types, counts_by_mode) for x in [feature]*matrix.shape[-1]
    ]

    feature_df = pd.DataFrame(
        {'feature_id' : np.arange(all_counts.shape[-1]).astype(str),
        'feature_type' : feature_type}
    )

    add_expression_to_dynframe(output_path, output_path, feature_df, 
        counts = all_counts)


def add_arguments(parser):

    parser.add_argument('--scaffold', '-i', type = str, required = True)
    parser.add_argument('--datasets', '-d', type = str, nargs = "+",
        required = True)
    parser.add_argument('--rd-means', '-mu', type = float, nargs = "+", required = True)
    parser.add_argument('--rd-stds', '-std', type = float, nargs = "+", required = True)
    parser.add_argument('--pure-states', '-s', type = str, nargs = "+", required = True)
    parser.add_argument('--feature-types', '-f', type = str, nargs = "+", required = True)
    parser.add_argument('--max-std', type = float, default = 2.5)
    parser.add_argument('--cell-state-col', '-col', type = str, default = 'leiden')
    parser.add_argument('--counts-layer', '-l', type = str, default = 'counts')
    parser.add_argument('--n-jobs', '-j', type = int, default = 1)
    parser.add_argument('--seed', type = int, default = None)
    parser.add_argument('--output-path', '-o', type = str, required = True)

