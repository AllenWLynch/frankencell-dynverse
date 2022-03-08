import mira
import numpy as np
from frankencell.preprocessing.basic_preprocessing import basic_rna_preprocessing, basic_atac_preprocessing
from frankencell.utils import read_dynframe, select_features, write_cell_info, add_dimred_prior, add_expression_to_dynframe
import logging
from shutil import copyfile

def init_rna_model(seed, kl_strategy, hidden = 128):
    rna_model = mira.topics.ExpressionTopicModel(
        beta=0.90,
        batch_size=32,
        seed = seed,
        encoder_dropout=0.05,
        num_topics = 8,
        decoder_dropout = 0.2,
        num_epochs = 60,
        counts_layer= 'counts',
        kl_strategy = kl_strategy,
	hidden = hidden
    )
    rna_model.set_learning_rates(1e-2, 2e-1)

    return rna_model
    

def init_atac_model(seed, kl_strategy, hidden = 128):

    atac_model = mira.topics.AccessibilityTopicModel(
        beta=0.93,
        batch_size=64,
        seed = seed,
        encoder_dropout=0.07,
        num_topics = 7,
        counts_layer= 'counts',
	hidden = hidden,
        kl_strategy = kl_strategy,
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
    training_args,
    seed = None,
    tune = True,
    min_cells = 25,
    min_dispersion = 0.7,
    kl_strategy = 'monotonic',
    hidden = 128
):  
    use_atac_features = not atac_data is None
    use_rna_features = not rna_data is None

    if use_rna_features:

        rna_model = init_rna_model(seed, kl_strategy, hidden = hidden)
        rna_data = basic_rna_preprocessing(rna_data, min_cells, min_dispersion)

        if tune:
            rna_model = tune_model(rna_data, rna_model, dataset_id+'_rna_study.pkl', **training_args)
        else:
            rna_model.fit(rna_data)

        rna_model.predict(rna_data)
        rna_model.save(dataset_id + '.rna-model.pth')
    
    if use_atac_features:

        basic_atac_preprocessing(atac_data, min_cells)
        
        atac_model = init_atac_model(seed, kl_strategy, hidden = hidden)

        if tune:
            atac_model = tune_model(atac_data, atac_model, dataset_id+'_atac_study.pkl', **training_args)
        else:
            atac_model.fit(atac_data)
            
        atac_model.predict(atac_data)
        atac_model.save(dataset_id + '.atac-model.pth')


    if use_rna_features and not use_atac_features:
        X = [rna_data.obsm['X_topic_compositions'], np.zeros((len(rna_data), 1))]
    elif use_atac_features and not use_rna_features:
        X = [np.zeros((len(rna_data), 1)), atac_data.obsm['X_topic_compositions']]
    else:
        X = [rna_data.obsm['X_topic_compositions'], np.zeros((len(rna_data), 1)), atac_data.obsm['X_topic_compositions']]
    rna_data.obsm['X_topic_compositions'] = np.hstack(X)
    return rna_data


def redundant_load(adata, feature_type, enforce_typing = True):
    try:
        data = select_features(adata, feature_type)
    except (IndexError, ValueError) as err:
        if enforce_typing:
            raise ValueError('Dataset not labeled with "feature_types" column, cannot separate RNA and ATAC features.')
        else:
            if isinstance(err, ValueError):
                logging.warning('Feature type {} not found in dataset, assuming all features are of type {}.'.format(feature_type, feature_type))
            else:
                logging.warning('Features not labeled with type, assuming all features are of type {}.'.format(feature_type, feature_type))

    return data
                
def main(
    dynframe_path,
    output_path = None,
    rna_feature_type = 'RNA',
    atac_feature_type = 'ATAC',
    use_atac_features = True,
    use_rna_features = True,
    tuning_iters = 32,
    min_topics = 7,
    max_topics = 13,
    max_dropout = 0.15,
    min_epochs = 30,
    max_epochs = 60,
    cv = 5,
    train_size=0.8,
    seed = None,
    tune = True,
    kl_strategy = 'monotonic',
    min_cells = 25,
    min_dispersion = 0.7,
    hidden = 128,
):

    training_args = dict(
        iters = tuning_iters, cv = cv, 
        min_epochs = min_epochs, max_epochs = max_epochs,
        min_topics = min_topics, max_topics = max_topics, max_dropout = max_dropout,
        train_size = train_size, 
    )

    adata = read_dynframe(dynframe_path)

    if use_rna_features:
        rna_data = redundant_load(adata, rna_feature_type, enforce_typing = use_atac_features)
    else:
        rna_data = None

    if use_atac_features:
        atac_data = redundant_load(adata, atac_feature_type, enforce_typing = use_rna_features)
    else:
        atac_data = None

    assert(not atac_data is None or not rna_data is None)

    out_data = train(
        rna_data, 
        atac_data,
        output_path,
        training_args,
        seed = seed,
        tune = tune,
        kl_strategy = kl_strategy,
        min_cells = min_cells,
        min_dispersion = min_dispersion,
	hidden = hidden
    )

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        out_data.var.reset_index(), 
        counts = out_data.layers['counts'], 
        expression = out_data.X.copy()
    )
    
    add_dimred_prior(output_path, out_data.obsm['X_topic_compositions'])
