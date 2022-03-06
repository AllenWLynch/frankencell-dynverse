import mira
import numpy as np
from frankencell.preprocessing.basic_preprocessing import basic_rna_preprocessing, basic_atac_preprocessing
from frankencell.utils import read_dynframe, select_features, write_cell_info, add_dimred_prior, add_expression_to_dynframe
import logging
from shutil import copyfile

def init_rna_model(seed):
    rna_model = mira.topics.ExpressionTopicModel(
        beta=0.90,
        batch_size=32,
        seed = seed,
        encoder_dropout=0.05,
        num_topics = 8,
        decoder_dropout = 0.2,
        num_epochs = 60,
        counts_layer= 'counts'
    )
    rna_model.set_learning_rates(1e-2, 2e-1)

    return rna_model
    

def init_atac_model(seed):

    atac_model = mira.topics.AccessibilityTopicModel(
        beta=0.93,
        batch_size=64,
        seed = seed,
        encoder_dropout=0.07,
        num_topics = 5,
        counts_layer= 'counts'
    )
    atac_model.set_learning_rates(1e-2, 2e-1)
    
    return atac_model


def tune_model(adata, model, save_name, **kwargs):

    train_size = kwargs.pop('train_size')

    tuner = mira.topics.TopicModelTuner(
        model, save_name = save_name,
        **kwargs
    )

    tuner.train_test_split(adata, train_size= train_size)
    tuner.tune(adata)

    best_model = tuner.select_best_model(adata)

    return best_model

def train(
    rna_data, 
    atac_data,
    dataset_id,
    seed = None,
    tune = True,
    box_cox = 0.5,
    min_cells = 25,
    min_dispersion = 0.7,
    **training_args,
):  
    use_atac_features = not atac_data is None

    rna_model = init_rna_model(seed)
    rna_data = basic_rna_preprocessing(rna_data, min_cells, min_dispersion)

    if tune:
        rna_model = tune_model(rna_data, rna_model, dataset_id+'_rna_study.pkl', **training_args)
    else:
        rna_model.fit(rna_data)

    rna_model.predict(rna_data)
    rna_model.get_umap_features(rna_data, box_cox = box_cox)
    
    if use_atac_features:

        basic_atac_preprocessing(atac_data, min_cells)
        
        atac_model = init_atac_model(seed)

        if tune:
            atac_model = tune_model(atac_data, atac_model, dataset_id+'_atac_study.pkl', **training_args)
        else:
            atac_model.fit(atac_data)
            
        atac_model.predict(atac_data)
        atac_model.get_umap_features(atac_data, box_cox = box_cox)

        rna_data, atac_data = mira.utils.make_joint_representation(
            rna_data, atac_data
        )

        rna_data.obsm['ATAC_topic_compositions'] = atac_data.obsm['X_topic_compositions']
        rna_data.obsm['ATAC_umap_features'] = atac_data.obsm['X_umap_features']

    return rna_data

        
def main(
    dynframe_path,
    output_path = None,
    rna_feature_type = 'RNA',
    atac_feature_type = 'ATAC',
    use_atac_features = True,
    tuning_iters = 32,
    min_topics = 5,
    max_topics = 12,
    max_dropout = 0.05,
    min_epochs = 20,
    max_epochs = 40,
    cv = 5,
    train_size=0.8,
    seed = None,
    tune = True,
    box_cox = 0.5,
    min_cells = 25,
    min_dispersion = 0.7
):

    training_args = dict(
        tuning_iters = tuning_iters, cv = cv, 
        min_epochs = min_epochs, max_epochs = max_epochs,
        min_topics = min_topics, max_topics = max_topics, max_dropout = max_dropout,
        train_size = train_size, 
    )

    adata = read_dynframe(dynframe_path)
    cell_info = {k : np.array(v) for k,v in adata.obs.reset_index().to_dict(orient = 'list').items()}

    try:
        rna_data = select_features(adata, rna_feature_type)
    except (IndexError, ValueError) as err:
        if use_atac_features:
            raise ValueError('Dataset not labeled with "feature_types" column, cannot separate RNA and ATAC features.')

        if isinstance(err, ValueError):
            logging.warning('Feature type {} not found in dataset, assuming all features are of type {}.'.format(rna_feature_type, rna_feature_type))
        else:
            logging.warning('Features not labeled with type, assuming all features are of type {}.'.format(rna_feature_type, rna_feature_type))

    if use_atac_features:
        atac_data = select_features(adata, atac_feature_type)
    else:
        atac_data = None

    out_data = train(
        rna_data, 
        atac_data,
        dynframe_path,
        seed = seed,
        tune = tune,
        box_cox = box_cox,
        min_cells = min_cells,
        min_dispersion = min_dispersion,
        **training_args,
    )

    save_cols = ['X_topic_compositions', 'X_umap_features']
    if use_atac_features:
        dimred = out_data.obsm['X_joint_umap_features']
        save_cols+=['ATAC_topic_compositions', 'ATAC_umap_features']
    else:
        dimred = out_data.obsm['X_umap_features']

    cell_info = {
        **cell_info,
        **{
            f + '_' + str(i) : out_data.obsm[f][:, i]
            for f in save_cols for i in range(out_data.obsm[f].shape[1])
        }
    }

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        out_data.var.reset_index(), 
        counts = out_data.layers['counts'], 
        expression = out_data.X.copy()
    )
    
    add_dimred_prior(output_path, dimred)

    write_cell_info(output_path, cell_info)