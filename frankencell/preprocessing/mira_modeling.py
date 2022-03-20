import mira
import numpy as np
from frankencell.preprocessing.basic_preprocessing import basic_rna_preprocessing, basic_atac_preprocessing
from frankencell.utils import read_dynframe, select_features, write_cell_info, add_dimred_prior, add_expression_to_dynframe
import logging
from shutil import copyfile
from sklearn.model_selection import ShuffleSplit
import os

def init_rna_model(model_args):

    default_args = dict(
        beta=0.92,
        batch_size=64,
        encoder_dropout=0.05,
        num_topics = 7,
        decoder_dropout = 0.2,
        num_epochs = 45,
    )

    default_args.update(model_args)

    rna_model = mira.topics.ExpressionTopicModel(
        counts_layer= 'counts',
        **default_args
    )
    rna_model.set_learning_rates(1e-2, 2e-1)

    return rna_model
    

def init_atac_model(model_args):

    default_args = dict(
        beta=0.92,
        batch_size=64,
        encoder_dropout=0.05,
        num_topics = 7,
        decoder_dropout = 0.2,
        num_epochs = 45,
    )

    default_args.update(model_args)

    atac_model = mira.topics.AccessibilityTopicModel(
        counts_layer= 'counts',
        **default_args
    )
    atac_model.set_learning_rates(5e-3, 1e-1)
    
    return atac_model


def tune_model(adata, model, save_name, tuning_args):

    default_args = dict(
        iters = 32,
        min_topics = 5,
        max_topics = 13,
        max_dropout = 0.1,
        min_epochs = 40,
        max_epochs = 60,
        cv = 'shufflesplit',
        train_size=0.9,
        top_n_trials = 5,
    )

    default_args.update(tuning_args)

    train_size = default_args.pop('train_size')
    cv = default_args.pop('cv')
    top_n_trials = default_args.pop('top_n_trials')

    if cv == 'shufflesplit':
        cv = ShuffleSplit(n_splits= 5, train_size=train_size)
    if cv == 'cv':
        cv = 5
    
    tuner = mira.topics.TopicModelTuner(
        model, save_name = save_name, cv = cv,
        **default_args
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

def parse_dataset_name(dataset_id, model_prefix, type):

    return os.path.join(os.path.dirname(dataset_id),
        '.'.join(os.path.basename(dataset_id).split('.')[0:2]) + '.' + model_prefix + '.h5.{}-model.pth'.format(type.lower()))

def train(
    rna_data, 
    atac_data,
    dataset_id,
    tuning_args,
    rna_model_args,
    atac_model_args,
    train_style = 'tune',
    min_cells = 25,
    min_dispersion = 0.7,
    rna_prefix = None,
    atac_prefix = None
):  
    use_atac_features = not atac_data is None
    use_rna_features = not rna_data is None

    if use_rna_features:

        rna_model = init_rna_model(rna_model_args)
        rna_data = basic_rna_preprocessing(rna_data, min_cells, min_dispersion)

        if train_style == 'tune':
            rna_model = tune_model(rna_data, rna_model, dataset_id+'_rna_study.pkl', tuning_args)
        elif train_style == 'fit':
            rna_model.fit(rna_data)
        elif train_style == 'load': 
            model_path = parse_dataset_name(dataset_id, rna_prefix, 'rna')
            rna_model = mira.topics.ExpressionTopicModel.load(model_path)
        else:
            raise ValueError('train_style="{}" not available.'.format(train_style))

        rna_model.predict(rna_data)
        rna_model.save(dataset_id + '.rna-model.pth')
    
    if use_atac_features:

        basic_atac_preprocessing(atac_data, min_cells)
        
        atac_model = init_atac_model(atac_model_args)

        if train_style == 'tune':
            atac_model = tune_model(atac_data, atac_model, dataset_id+'_atac_study.pkl', tuning_args)
        elif train_style == 'fit':
            atac_model.fit(atac_data)
        elif train_style == 'load': 
            model_path = parse_dataset_name(dataset_id, atac_prefix, 'atac')
            atac_model = mira.topics.AccessibilityTopicModel.load(model_path)
        else:
            raise ValueError('train_style="{}" not available.'.format(train_style))
            
        atac_model.predict(atac_data)
        atac_model.save(dataset_id + '.atac-model.pth')


    if use_rna_features and not use_atac_features:
        X = [rna_data.obsm['X_topic_compositions'], np.zeros((len(rna_data), 1))]
        rna_data.obsm['X_topic_compositions'] = np.hstack(X)
        return rna_data
    elif use_atac_features and not use_rna_features:
        X = [np.zeros((len(atac_data), 1)), atac_data.obsm['X_topic_compositions']]
        atac_data.obsm['X_topic_compositions'] = np.hstack(X)
        return atac_data
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
    train_style = 'tune',
    min_cells = 25,
    min_dispersion = 0.7,
    tuning_args = dict(),
    rna_model_args = dict(),
    atac_model_args = dict(),
    rna_prefix = None,
    atac_prefix = None
):

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
        tuning_args,
        rna_model_args,
        atac_model_args,
        train_style = train_style,
        min_cells = min_cells,
        min_dispersion = min_dispersion,
        rna_prefix = rna_prefix,
        atac_prefix = atac_prefix,
    )

    add_expression_to_dynframe(
        dynframe_path, 
        output_path,
        out_data.var.reset_index(), 
        counts = out_data.layers['counts'], 
        expression = out_data.X.copy()
    )
    
    add_dimred_prior(output_path, out_data.obsm['X_topic_compositions'])
