from random import choices
from dynclipy.read import ro as dynro
import os

def main(*, goldstandard, test_dataset, run_file, definition_file,
    method_output_path, results_output_path,
    param_string = ''):

    path = os.path.dirname(os.path.abspath(__file__))
    rstring = '''
goldstandard <- dynutils::read_h5("{goldstandard}")
model <- dynutils::read_h5("{test_dataset}")

test_method <- dynwrap::create_ti_method_definition(
    "{definition_file}",
    "{run_file}",
)

model <- dynwrap::infer_trajectory(model, test_method({param_string}), verbose = TRUE,
            give_priors = c('start_id','end_id','dimred'))

model <- dynwrap::add_cell_waypoints(model)

results <- dyneval::calculate_metrics(goldstandard, model,
                        metrics = c("correlation","F1_branches","edge_flip"))

dynutils::write_h5(model, "{method_output_path}")

write.table(t(results), file = "{results_output_path}", 
    sep = "\t", col.names = FALSE, quote = FALSE)

    '''.format(
        path = path,
        run_file = run_file,
        definition_file = definition_file,
        goldstandard = goldstandard,
        test_dataset = test_dataset,
        method_output_path = method_output_path,
        param_string = param_string,
        results_output_path = results_output_path
    )

    dynro.r(rstring)

def add_arguments(parser):
    parser.add_argument('--truth-dataset', '-t', type = str, required = True)
    parser.add_argument('--test-dataset', type = str, required = True)
    parser.add_argument('--method-outfile', '-o', type = str, required = True)
    parser.add_argument('--results-outfile', '-r', type = str, required = True)
    parser.add_argument('--run-file', type = str, required = True)
    parser.add_argument('--definition-file', type = str, required = True)
    parser.add_argument('--parameters', '-p', type = str, required = False,
        default = '')