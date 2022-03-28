
import argparse
import sys
import yaml
from yaml import Loader
from .generate_cells import generate_frankentrajectory
from .generate_cells import add_arguments as scaffold_add_arguments

from .mix_cells import mix_frankencells, read_datasets
from .mix_cells import add_arguments as mixer_add_arguments

from .preprocessing.pca import main as pca_preprocessing
from .preprocessing.pca import add_arguments as pca_add_arguments

from .preprocessing.lsi import main as lsi_preprocessing
from .preprocessing.lsi import add_arguments as lsi_add_arguments

from .evaluate_method import main as evaluate
from .evaluate_method import add_arguments as eval_add_arguments

from .get_metrics import main as get_metrics
from .get_metrics import add_arguments as metrics_add_arguments

from .main import main as maketests
from .main import add_arguments as maketests_add_arguments

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
subparsers = parser.add_subparsers(help = 'commands')

def maketests_wrapper(args):

    config = yaml.load(args.config, Loader)

    maketests(
        config = config,
        restart= args.restart,
        mem = args.mem_mb,
        snake_args=args.snake_args
    )

maketests_subparser = subparsers.add_parser('gen')
maketests_add_arguments(maketests_subparser)
maketests_subparser.set_defaults(func = maketests_wrapper)

'''def generate_scaffold(args):

    generate_frankentrajectory(
        args.state_composition,
        args.outfile,
        branch_times = args.branch_times,
        n_cells = args.n_cells,
        gamma = args.gamma,
        max_depth = args.max_depth,
        max_width = args.max_width,
        ptime_alpha = args.ptime_alpha,
        ptime_beta = args.ptime_beta,
        sigmoid_approach = not args.no_sigmoid,
        sigmoid_aggression = args.sigmoid_aggression,
        seed = args.seed,
    )

scaffold_subparser = subparsers.add_parser('scaffold', fromfile_prefix_chars='@')
scaffold_add_arguments(scaffold_subparser)
scaffold_subparser.set_defaults(func = generate_scaffold)

def mix_cells(args):

    mix_frankencells(
        dynframe = args.scaffold,
        datasets = read_datasets(args.datasets), 
        rd_means = args.rd_means, 
        rd_stds = args.rd_stds, 
        pure_states = args.pure_states, 
        feature_types = args.feature_types,
        max_std = args.max_std, 
        cell_state_col = args.state_col,
        counts_layer = args.counts_layer,
        n_jobs = args.n_jobs,
        seed = args.seed,
        output_path = args.outfile, 
    )

mix_cells_subparser = subparsers.add_parser('mix-cells', fromfile_prefix_chars='@')
mix_cells_subparser.fromfile_prefix_chars = '@'
mixer_add_arguments(mix_cells_subparser)
mix_cells_subparser.set_defaults(func = mix_cells)

def run_pca_preprocess(args):

    pca_preprocessing(
        args.dynframe,
        output_path=args.outfile,
        feature_type=args.feature_type,
        min_cells=args.min_cells,
        min_dispersion=args.min_dispersion
    )

pca_parser = subparsers.add_parser('pca-preprocess')
pca_add_arguments(pca_parser)
pca_parser.set_defaults(func = run_pca_preprocess)


def run_lsi_preprocess(args):

    lsi_preprocessing(
        args.dynframe,
        output_path=args.outfile,
        feature_type=args.feature_type,
        min_cells=args.min_cells,
    )

lsi_parser = subparsers.add_parser('lsi-preprocess')
lsi_add_arguments(lsi_parser)
lsi_parser.set_defaults(func = run_lsi_preprocess)

def run_evaluate(args):

    evaluate(
        goldstandard=args.truth_dataset,
        test_dataset = args.test_dataset,
        run_file = args.run_file,
        definition_file=  args.definition_file,
        method_output_path = args.method_outfile,
        results_output_path= args.results_outfile,
        param_string=args.parameters,
    )

evaluate_parser = subparsers.add_parser('evaluate')
eval_add_arguments(evaluate_parser)
evaluate_parser.set_defaults(func = run_evaluate)


def run_get_metrics(args):

    get_metrics(
        goldstandard=args.truth_dataset,
        test_dataset = args.test_dataset,
        results_output_path= args.results_outfile,
    )

metrics_subparser = subparsers.add_parser('get-metrics')
metrics_add_arguments(metrics_subparser)
metrics_subparser.set_defaults(func = run_get_metrics)'''


def main():
    #____ Execute commands ___

    args = parser.parse_args()

    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        args.func(args)