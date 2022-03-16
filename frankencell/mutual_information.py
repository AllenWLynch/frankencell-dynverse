
from asyncore import read
import numpy as np
from .utils import read_cell_info


def get_topics(data):

    topics = data.obsm['dimred']
    divider_col = np.argwhere(np.isclose(topics, 0).all(0))[0,0]

    if use_rep == 'joint' or use_rep == 'joint_constrained':
    if divider_col == 0 or divider_col == topics.shape[-1]-1:
        raise ValueError('User wanted to use "joint" representation, but only one representation was given')

    rep1, rep2 = topics[:, :divider_col], topics[:,divider_col+1:]

        if use_rep == 'joint':
            data.obsm['embedding_features'] = np.hstack([
                get_umap_features(rep1, box_cox), 
                get_umap_features(rep2, box_cox)]
            )
        else:
            data.obsm['embedding_features'] = np.hstack([
                l2_norm(get_umap_features(rep1, box_cox)), 
                l2_norm(get_umap_features(rep2, box_cox))]
            )

    elif use_rep == 'RNA':
        data.obsm['embedding_features'] = get_umap_features(topics[:, :divider_col], box_cox)
    elif use_rep == 'ATAC':
        data.obsm['embedding_features'] = get_umap_features(topics[:, divider_col+1:], box_cox)
    else:
        raise ValueError('Representation {} is unknown.'.format(use_rep))

    #print(data.obsm['embedding_features'], data.obsm['embedding_features'].shape, data.obs_names)

else:
    data.obsm['embedding_features'] = data.obsm['dimred']


def _mutual_information(x,y):

    x_marg = x.mean(0,keepdims = True)
    y_marg = y.mean(0, keepdims = True)

    joint = (x[:, np.newaxis, :] * y[:,:, np.newaxis])
    marg = (x_marg * y_marg.T)[np.newaxis, :,:]

    mutual_information = np.sum(joint*np.log2(joint/marg), axis = (-2,-1)).mean()

    return mutual_information


def get_topic_MI(data):

    cell_info = read_cell_info(data)

    print(cell_info.columns)
    
    topic_comps = cell_info.iloc[:, cell_info.columns.str.startswith('topic_')].values
    mixing_weights = cell_info.iloc[:, cell_info.columns.str.startswith('mix_weight_')].values

    print(topic_comps)
    print(mixing_weights)

    return _mutual_information(topic_comps, mixing_weights)

def main(*,data):

    print(get_topic_MI(data = data))


def get_arguments(parser):

    parser.add_arguments('--dynframe', '-d', type = str, required = True)