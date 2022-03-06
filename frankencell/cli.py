
import argparse
import sys
from .generate_cells import generate_frankentrajectory
from .generate_cells import add_arguments as scaffold_add_arguments

def generate_scaffold(args):

    generate_frankentrajectory(
        args.state_compositions,
        args.outfile,
        branch_times = args.branch_times,
        n_cells = args.n_cells,
        gamma = args.gamma,
        max_depth = args.max_depth,
        max_width = args.max_width,
        ptime_alpha = args.ptime_alpha,
        ptime_beta = args.ptime_beta,
        sigmoid_approach = args.sigmoid_approach,
        sigmoid_aggression = args.sigmoid_aggression,
        seed = args.seed,
    )

parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
subparsers = parser.add_subparsers(help = 'commands')

scaffold_subparser = subparsers.add_parser('scaffold')
scaffold_add_arguments(scaffold_subparser)
scaffold_subparser.set_defaults(func = generate_scaffold)

def main():
    #____ Execute commands ___
    args = parser.parse_args()

    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        args.func(args)