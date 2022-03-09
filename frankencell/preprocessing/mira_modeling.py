import mira
import numpy as np
from frankencell.preprocessing.basic_preprocessing import basic_rna_preprocessing, basic_atac_preprocessing
from frankencell.utils import read_dynframe, select_features, write_cell_info, add_dimred_prior, add_expression_to_dynframe
import logging
from shutil import copyfile
from sklearn.model_selection import ShuffleSplit

def init_rna_model(model_args):
    rna_model = mira.topics.ExpressionTopicModel(
        beta=0.90,
        batch_size=32,
        encoder_dropout=0.05,
        num_topics = 8,
        decoder_dropout = 0.2,
        num_epochs = 60,
        counts_layer= 'counts',
        **model_args
    )
    rna_model.set_learning_rates(1e-2, 2e-1)

    return rna_model
    

def init_atac_model(model_args):

    atac_model = mira.topics.AccessibilityTopicModel(
        beta=0.93,
        batch_size=32,
        encoder_dropout=0.07,
        num_topics = 7,
        counts_layer= 'counts',
        **model_args
    )
    atac_model.set_learning_rates(1e-2, 2e-1)
    
    return atac_model


def tune_model(adata, model, save_name, **kwargs):

    train_size = kwargs.pop('train_size')
    cv = kwargs.pop('cv')
    top_n_trials = kwargs.pop('top_n_trials')

    if cv == 'shufflesplit':
        cv = ShuffleSplit(n_splits= 5, train_size=train_size)
    if cv == 'cv':
        cv = 5

    tuner = mira.topics.TopicModelTuner(
        model, save_name = save_name, cv = cv,
        **kwargs
    )

    if top_n_trials == 1:
        adata.obs[tuner.test_column] = False
        adata.obs[tuner.test_column][np.random.choice(len(adata), size = 2)] = True
    else:
        tuner.train_test_split(adata, train_size= train_size)

    tuner.tune(adata)

    if top_n_trials == 1:
        best_model = model.set_params(**tuner.get_best_params(1)[0]).fit(adata)
    else:
        best_model = tuner.select_best_model(adata, top_n_trials=top_n_trials)

    return best_model

def train(
    rna_data, 
    atac_data,
    dataset_id,
    training_args,
    model_args,
    train_style = 'tune',
    min_cells = 25,
    min_dispersion = 0.7,
):  
    use_atac_features = not atac_data is None
    use_rna_features = not rna_data is None

    if use_rna_features:

        rna_model = init_rna_model(model_args)
        rna_data = basic_rna_preprocessing(rna_data, min_cells, min_dispersion)

        if train_style == 'tune':
            rna_model = tune_model(rna_data, rna_model, dataset_id+'_rna_study.pkl', **training_args)
        elif train_style == 'fit':
            rna_model.fit(rna_data)
        elif train_style == 'load': 
            rna_model = mira.topics.ExpressionTopicModel.load(dataset_id + '.rna-model.pth')
        else:
            raise ValueError('train_style="{}" not available.'.format(train_style))

        rna_model.predict(rna_data)
        rna_model.save(dataset_id + '.rna-model.pth')
    
    if use_atac_features:

        basic_atac_preprocessing(atac_data, min_cells)
        
        atac_model = init_atac_model(model_args)

        if train_style == 'tune':
            atac_model = tune_model(atac_data, atac_model, dataset_id+'_atac_study.pkl', **training_args)
        elif train_style == 'fit':
            atac_model.fit(atac_data)
        elif train_style == 'load': 
            atac_model = mira.topics.AccessibilityTopicModel.load(dataset_id + '.atac-model.pth')
        else:
            raise ValueError('train_style="{}" not available.'.format(train_style))
            
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
    train_style = 'tune',
    min_cells = 25,
    min_dispersion = 0.7,
    kl_strategy = 'monotonic',
    hidden = 128,
    layers = 3,
    top_n_trials = 5,
):

    training_args = dict(
        iters = tuning_iters, cv = cv, 
        min_epochs = min_epochs, max_epochs = max_epochs,
        min_topics = min_topics, max_topics = max_topics, max_dropout = max_dropout,
        train_size = train_size, top_n_trials = top_n_trials,
    )

    model_args = dict(
        kl_strategy = kl_strategy, hidden = hidden, seed = seed,
        num_layers = layers, 

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
        model_args,
        train_style = train_style,
        min_cells = min_cells,
        min_dispersion = min_dispersion,
    )

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        out_data.var.reset_index(), 
        counts = out_data.layers['counts'], 
        expression = out_data.X.copy()
    )
    
    add_dimred_prior(output_path, out_data.obsm['X_topic_compositions'])
