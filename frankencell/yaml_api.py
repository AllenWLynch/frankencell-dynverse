import argparse
from .generate_cells import generate_frankentrajectory
from .mix_cells import mix_frankencells, read_datasets
from .preprocessing.pca import main as pca_preprocessing
from .preprocessing.lsi import main as lsi_preprocessing
from .preprocessing.mira_modeling import main as mira_preprocessing
from .evaluate_method import main as evaluate
import yaml
from yaml import Loader

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
mix_cells_subparser.set_defaults(func = mix_frankencells)

pca_parser = subparsers.add_parser('pca-preprocess')
add_arguments(pca_parser)
pca_parser.set_defaults(func = pca_preprocessing)

lsi_parser = subparsers.add_parser('lsi-preprocess')
add_arguments(lsi_parser)
lsi_parser.set_defaults(func = lsi_preprocessing)

mira_parser = subparsers.add_parser('mira-preprocess')
add_arguments(mira_parser)
mira_parser.set_defaults(func = mira_preprocessing)

evaluate_parser = subparsers.add_parser('evaluate')
add_arguments(evaluate_parser)
evaluate_parser.set_defaults(func = evaluate)


def main():

    args = parser.parse_args()


    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        with open(args.yaml, 'r') as f:
            params = yaml.load(f, Loader)
            args.func(**params)
