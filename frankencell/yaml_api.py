import argparse
from .generate_cells import generate_frankentrajectory
from .mix_cells import mix_frankencells, read_datasets
import yaml
from yaml import Loader
import sys

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help = 'commands')

def add_arguments(parser):
    parser.add_argument('yaml', type = str)
    #parser.add_argument('--threads', default = 1, required = False, type = int)

scaffold_subparser = subparsers.add_parser('scaffold')
add_arguments(scaffold_subparser)
scaffold_subparser.set_defaults(func = generate_frankentrajectory)

mix_cells_subparser = subparsers.add_parser('mix-cells')
add_arguments(mix_cells_subparser)
mix_cells_subparser.add_argument('--n-jobs', '-j', type = int, default=1)
mix_cells_subparser.set_defaults(func = mix_frankencells)


def main():

    args = parser.parse_args()


    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        with open(args.yaml, 'r') as f:
            params = yaml.load(f, Loader)
        
        try:
            params['n_jobs'] = args.n_jobs
        except AttributeError:
            pass

        args.func(**params)
