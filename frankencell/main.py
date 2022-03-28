from typing_extensions import Required
import os
import subprocess
import argparse
import yaml

def write_scaffold_config(config, prefix, scaffold, SCAFFOLD_CONFIG, SCAFFOLD_PATH):
        
    with open(SCAFFOLD_CONFIG.format(
                prefix=prefix, scaffold=scaffold),'w') as f:
        yaml.dump({
            **config['default_scaffold_parameters'],
            **config['scaffolds'][scaffold],
            'output_path' : SCAFFOLD_PATH.format(scaffold=scaffold)
        }, f)


def write_test_config(config, prefix, scaffold, test, replicate, SCAFFOLD_PATH, TEST_CONFIG, TEST_PATH):
    
    outfile = TEST_PATH.format(scaffold=scaffold, test=test, rep=replicate)

    with open(TEST_CONFIG.format(
                prefix=prefix, scaffold=scaffold, test=test, rep=replicate),'w') as f:

        yaml.dump({
            **config['default_test_parameters'],
            **config['test_conditions'][test],
            'scaffold' : SCAFFOLD_PATH.format(scaffold=scaffold),
            'output_path' : outfile,
            'seed' : replicate + 1775,
        }, f)

    return outfile


def set_up_trials(config):

    prefix = config['prefix']
    SCAFFOLD_PATH = prefix+'scaffolds/{scaffold}.h5'
    TEST_PATH = prefix+'tests/{scaffold}.{test}-{rep}.h5'

    configdir_path = prefix + 'configs/'
    if not os.path.isdir(configdir_path):
        os.mkdir(configdir_path)
        
    SCAFFOLD_CONFIG = configdir_path + '{scaffold}.yaml'
    TEST_CONFIG = configdir_path + '{scaffold}.{test}-{rep}.yaml'

    targets = []
    for scaffold in config['scaffolds']:
        write_scaffold_config(config, prefix, scaffold, SCAFFOLD_CONFIG, SCAFFOLD_PATH)
        for test in config['test_conditions']:
            for replicate in range(1, config['test_replicates'] + 1):
                targets.append(
                    write_test_config(config, prefix, scaffold, test, replicate, SCAFFOLD_PATH, TEST_CONFIG, TEST_PATH)
                )

    snake_config = {
        'SCAFFOLD_PATH' : SCAFFOLD_PATH,
        'SCAFFOLD_CONFIG' : SCAFFOLD_CONFIG,
        'TEST_PATH' : TEST_PATH,
        'TEST_CONFIG' : TEST_CONFIG,
    }

    return targets, snake_config


def run_pipeline(*, snake_config, mem, snake_args, targets):

    snakefile_path = os.path.join(
        os.path.dirname(__file__), 'pipeline', 'snakefile_gen'
    )

    snake_args = [
        'snakemake',
        *targets,
        '-s', snakefile_path,
        '--resources', 'local_mem={}'.format(str(mem)),
        '--config', *['{}="{}"'.format(k, v) for k,v in snake_config.items()],
        *snake_args,
    ]

    subprocess.run(
        snake_args
    )


def add_arguments(parser):
    parser.add_argument('config', type = argparse.FileType('r'), help = 'Test configuration file.')
    parser.add_argument('--mem_mb', '-mem', required = True, type = int, help = 'Local memory available on your system. Changes how many concurrent jobs may run.')
    parser.add_argument('--restart', '-r', default = False, type = bool, help = 'Rerun the entire pipeline from the start.')
    parser.add_argument('--snake_args', '-s', nargs = argparse.REMAINDER, help = 'arguments to pass to snakemake')

def main(*, config, restart, mem, snake_args):
    
    targets, snake_config = set_up_trials(config)
    snake_config['restart'] = restart

    run_pipeline(
        targets = targets,
        snake_config = snake_config, 
        mem = mem, 
        snake_args = snake_args
    )