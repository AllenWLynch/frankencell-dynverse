import yaml

with open('frankencell/pipeline/config.yaml', 'r') as f:
    config = yaml.load(f, yaml.Loader)

def write_scaffold_config(config, prefix, scaffold):
    
    with open('{prefix}.{scaffold}.yaml'.format(
                prefix=prefix, scaffold=scaffold),'w') as f:
        yaml.dump({
            **config['default_scaffold_parameters'],
            **config['scaffolds'][scaffold],
            'output_path' : '{prefix}.{scaffold}.h5'.format(prefix=prefix, scaffold=scaffold)
        }, f)


def write_test_config(config, prefix, scaffold, test, replicate):
    
    with open('{prefix}.{scaffold}.{test}-{rep}.yaml'.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate),'w') as f:
        yaml.dump({
            **config['default_test_parameters'],
            **config['test_conditions'][test],
            'scaffold' : '{prefix}.{scaffold}.h5'.format(prefix = prefix, scaffold = scaffold),
            'output_path' : '{prefix}.{scaffold}.{test}-{rep}.h5'.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate),
            'seed' : replicate + 1775,
        }, f)
        
def write_preprocess_config(config, prefix, 
                        scaffold, test, replicate, preprocess):

    with open('{prefix}.{scaffold}.{test}-{rep}.{preprocess}.yaml'.format(
                prefix=prefix, scaffold=scaffold, test=test, 
                preprocess=preprocess, rep=replicate),'w') as f:
        yaml.dump({
            **config['default_preprocessing_parameters'],
            **config['preprocessing'][preprocess]['parameters'],
            'dynframe_path' : '{prefix}.{scaffold}.{test}-{rep}.h5'.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate),
            'output_path' : '{prefix}.{scaffold}.{test}-{rep}.{preprocess}.h5'.format(
                prefix=prefix, scaffold=scaffold, 
                test=test, preprocess=preprocess, rep=replicate)
        }, f)

def write_method_config(config, prefix, scaffold, test, replicate, method,
                       trial):
    
    preprocess = config['methods'][method]['trials'][trial]['preprocessing']

    write_scaffold_config(config, prefix, scaffold)
    write_test_config(config, prefix, scaffold, test, replicate)
    write_preprocess_config(config, prefix, scaffold, test, replicate,
                           preprocess)
    
    result_file = '{prefix}.{scaffold}.{test}-{rep}.{preprocess}.{method}.{trial}.results.csv'.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate,
                preprocess=preprocess, method = method, trial = trial)
    
    with open('{prefix}.{scaffold}.{test}-{rep}.{preprocess}.{method}.{trial}.yaml'.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate,
                preprocess=preprocess, method=method, trial=trial),'w') as f:
        yaml.dump({**config['default_method_parameters'],
         'run_file' : config['methods'][method]['run_file'],
         'definition_file' : config['methods'][method]['definition_file'],
         'goldstandard' : '{prefix}.{scaffold}.h5'.format(prefix=prefix, scaffold=scaffold),
         'test_dataset' : '{prefix}.{scaffold}.{test}-{rep}.{preprocess}.h5'.format(
                prefix=prefix, scaffold=scaffold, 
                test=test, preprocess=preprocess, rep=replicate),
         'method_output_path' : '{prefix}.{scaffold}.{test}-{rep}.{preprocess}.{method}.{trial}.h5'.format(
                prefix=prefix, scaffold=scaffold, test=test, 
                preprocess=preprocess, method = method, 
                trial = trial, rep=replicate),
         'results_output_path' : result_file,
         'param_string' : config['methods'][method]['trials'][trial]['parameters']
        }, f)
        
    return result_file


def set_up_trials(config):

    results_files = []
    prefix = config['prefix']
    for scaffold in config['scaffolds']:
        for test in config['test_conditions']:
            for replicate in range(1, config['test_replicates'] + 1):
                for method in config['methods']:
                    for trial in config['methods'][method]['trials']:
                        #print(prefix, scaffold, test, replicate, method, trial)
                        results_files.append(
                            write_method_config(
                                config, prefix, scaffold, test, replicate, method, trial
                            )
                        )

    return results_files
