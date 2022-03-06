from dynclipy.read import ro as dynro
import os

def main(*, goldstandard, model, method, outfile,
    param_string = ''):

    path = os.path.dirname(os.path.abspath(__file__))
    rstring = '''
mira <- dynwrap::create_ti_method_definition(
    "{path}/methods/mira/definition.yaml",
    "{path}/methods/mira/run",
)

lsi <- dynwrap::create_ti_method_definition(
    "{path}/methods/signac/definition.yaml",
    "{path}/frankencell/methods/signac/run.R",
)

sling <- dynwrap::create_ti_method_definition(
    "{path}/methods/custom_slingshot/definition.yaml",
    "{path}/methods/custom_slingshot/run.R",
)

goldstandard <- dynutils::read_h5("{goldstandard}")
model <- dynutils::read_h5("{model}")

methods_list <- c(mira, sling, lsi)
names(methods_list) <- c("mira","slingshot", "lsi")

test_method <- methods_list[["{method}"]]

model <- dynwrap::infer_trajectory(model, test_method({param_string}), verbose = TRUE,
            give_priors = c('start_id','end_id','dimred'))

model <- dynwrap::add_cell_waypoints(model)

results <- dyneval::calculate_metrics(goldstandard, model,
                        metrics = c("correlation","F1_branches","edge_flip"))

write.csv(results, file = "{outfile}")
    '''.format(
        path = path,
        goldstandard = goldstandard,
        model = model,
        method = method,
        outfile = outfile,
        param_string = param_string
    )
    
    dynro.r(rstring)