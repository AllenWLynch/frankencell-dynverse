from dynclipy.read import ro as dynro
import os
from sklearn.metrics import fbeta_score
import pandas as pd


def get_modified_F1_score(F1_output_path):

    assignments = pd.read_csv(F1_output_path, sep = ' ').set_index('cell_id').rename(
        columns  = {'group_prediction' : 'prediction', 'group_dataset' : 'actual'}
    )

    mappings = assignments.groupby(['actual','prediction']).size().reset_index().sort_values(0).groupby('prediction').tail(1)

    mappings = dict(zip(mappings.actual, mappings.prediction))

    actual_branches = list(mappings.keys())

    leftout = set(assignments.actual).difference(set(actual_branches))

    if len(leftout) > 0:

        def get_startnode(x):
            return x.split('-')[0]

        def get_endnode(x):
            return x.split('>')[1]
        
        added_leg = False
        for branch in leftout:
            for start_node, in_branch in zip(map(get_startnode, actual_branches), actual_branches): 
                if get_endnode(branch) == start_node:
                    mappings[branch] = mappings[in_branch]
                    added_leg = True

        if not added_leg:
            for branch in leftout:
                potential_legs = []
                
                for end_node, in_branch in zip(map(get_endnode, actual_branches), actual_branches): 
                    if get_startnode(branch) == end_node:
                        potential_legs.append(branch)
                
                if len(potential_legs) > 0:
                    add_leg = assignments.groupby('actual').count().loc[potential_legs].sort_values('prediction').head(1).index.values[0]
                    mappings[branch] = mappings[add_leg]
                    added_leg = True
                    break

    mapped_actual = assignments.actual.map(mappings).fillna('None').values

    return fbeta_score(mapped_actual, assignments.prediction.values, average='macro', beta=1)



def main(*, goldstandard, test_dataset, run_file, definition_file,
    method_output_path, results_output_path,
    param_string = ''):

    path = os.path.dirname(os.path.abspath(__file__))
    F1_output_path = os.path.join(
        os.path.dirname(results_output_path),
        'F1_data.csv',
    )

    rstring = '''

library(tibble)

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

groups_dataset <- goldstandard %>% dynwrap::group_onto_trajectory_edges()
groups_prediction <- model %>% dynwrap::group_onto_trajectory_edges()
groups_dataset <- groups_dataset %>% as.character() %>% enframe("cell_id", "group_dataset")
groups_dataset_levels <- unique(groups_dataset$group_dataset) %>% na.omit()
groups_prediction <- groups_prediction %>% as.character() %>% enframe("cell_id", "group_prediction")
groups_prediction_levels <- unique(groups_prediction$group_prediction) %>% na.omit()
groups <- full_join(groups_dataset, groups_prediction, "cell_id")
write.table(groups, "{F1_output_path}", sep = "\t")

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
        results_output_path = results_output_path,
        F1_output_path = F1_output_path,
    )

    dynro.r(rstring)

    mod_f1_score = get_modified_F1_score(F1_output_path)

    with open(results_output_path, 'a') as f:
        print('mod_F1_branches', mod_f1_score, sep = '\t')

    os.remove(F1_output_path)
    

def add_arguments(parser):
    parser.add_argument('--truth-dataset', '-t', type = str, required = True)
    parser.add_argument('--test-dataset', type = str, required = True)
    parser.add_argument('--method-outfile', '-o', type = str, required = True)
    parser.add_argument('--results-outfile', '-r', type = str, required = True)
    parser.add_argument('--run-file', type = str, required = True)
    parser.add_argument('--definition-file', type = str, required = True)
    parser.add_argument('--parameters', '-p', type = str, required = False,
        default = '')