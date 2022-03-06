import numpy as np
import argparse
import pandas as pd
import dynclipy as dyn
import networkx as nx
from networkx.generators.classic import _tree_edges
import h5py as h5
import logging
from .utils import write_cell_info

from dynclipy.dataset import add_adder, Dataset
add_adder(Dataset, "add_prior_information")
add_adder(Dataset, "add_cell_waypoints")

class CRP:
    
    def __init__(self, gamma, max_depth, max_width = 2):
        self.gamma = gamma
        self.max_depth = max_depth
        self.max_width = max_width
        self.tables = []
        self.n = 0
    
    def _nested_sum(self, table):
        total = 0
        if not isinstance(table, list):
            return table
        for i in table:
            if isinstance(i, list):
                total += self._nested_sum(i)
            else:
                total += i
        return total
        
    def _seat(self, table, depth):
        
        is_leaf= depth+1 >= self.max_depth
               
        n_seated = np.array([self._nested_sum(t) for t in table])
        n_total = n_seated.sum() + 1
        
        if n_total == 0:
            selected_table = 0
        else:
            gamma = self.gamma
            if len(table) == self.max_width:
                gamma = 0
                
            denom = gamma + n_total -1
            p_table = np.concatenate([n_seated/denom, [gamma/denom]], axis = -1)
                        
            selected_table = np.random.choice(range(len(table)+1), p = p_table)
        
        recurse_table = []
        if selected_table == len(table):
                        
            if is_leaf:
                table.append(1)
                return [selected_table]
            else:
                table.append(recurse_table)
                return [selected_table, *self._seat(recurse_table, depth+1)]
                     
        else:
            recurse_table = table[selected_table]
        
            if is_leaf:
                return [selected_table]
            else:
                return [selected_table, *self._seat(recurse_table, depth+1)]
    
    def new_customer(self):
        return [0,*self._seat(self.tables, 1)]
    
    
def sigmoid(x):
    return 1/(1+np.exp(-x))


def get_idx_from_path(path):
    idx = [1]
    for p in path[1:]:
        idx.append(2*idx[-1] + p)
        
    return idx

def get_mixing_weights(*,state_compositions, paths, transition_mixing):

    paths = np.array([[0, *get_idx_from_path(p)] for p in paths])
    transition_mixing = transition_mixing[:,:,np.newaxis]

    return (state_compositions[paths.astype(int)] * transition_mixing).sum(1)

def generate_branching_process(
    branch_times,
    n_cells = 1000,
    gamma = 0.1,
    max_depth = 2,
    max_width = 2,
    ptime_alpha = 1.,
    ptime_beta = 1.,
    sigmoid_approach = True,
    sigmoid_aggression = 5):

    assert(len(branch_times) == (max_depth +1))
    branch_times = np.array(branch_times)
    crp = CRP(gamma, max_depth, max_width=max_width)
    
    min_sigmoid, max_sigmoid = sigmoid(-0.5 * sigmoid_aggression), sigmoid(0.5 * sigmoid_aggression)
    
    cells = []
    for _ in range(n_cells):
        
        pseudotime = np.random.beta(ptime_alpha, ptime_beta)
        path = crp.new_customer()
                        
        level = np.argmin(pseudotime > branch_times) - 1
        branch_len = branch_times[level+1] - branch_times[level]
        progress = (pseudotime - branch_times[level])/(branch_len)
        
        if sigmoid_approach:
            x = sigmoid_aggression*(progress - 0.5)
            progress = (sigmoid(x) - min_sigmoid)/(max_sigmoid - min_sigmoid)
            
        mixing = np.concatenate([np.zeros(level), 
                                 [1-progress, progress],
                                np.zeros(len(path) - level - 1)])
        
        state = get_idx_from_path(path[:level+1])[-1] - 1

        cells.append((pseudotime, path, mixing, state,
            progress * branch_len + branch_times[level], progress))

    random_order = np.random.permutation(len(cells))
    cells = [cells[j] for j in random_order]
        
    return list(
        map(np.array, list(zip(*cells))) #pseudotime, path, mixing, state, ?, progress
    )

def make_graph(states):

    n_nodes = len(np.unique(states)) + 1
    G = nx.empty_graph(n_nodes, create_using=nx.DiGraph)
    G.add_edges_from([(-1,0)])
    G.add_edges_from(_tree_edges(n_nodes - 1, 2))
    G.remove_node(n_nodes-1)
    
    return G

def simplify_tree(tree, state_compositions, states, progress):
    
    for node in nx.dfs_postorder_nodes(tree, -1):
        
        children = list(nx.dfs_postorder_nodes(tree, source = node))[:-1]
        child_states = [x+1 for x in children]
        num_children = len(children)
        
        if num_children > 1 and np.isclose(state_compositions[child_states], state_compositions[child_states[0]]).all():
            states = np.where(
                np.isin(states, children[1:]), children[0], states
            )
            
            tree.remove_nodes_from(children[1:])
            num_children-=len(children[1:])
            
        
        if num_children == 1:
            progress = np.where(states ==  children[0], 0.5*progress + 0.5, progress)
            progress = np.where(states == node, 0.5*progress , progress)
            states = np.where(states == children[0], node, states)
        
            tree.remove_node(children[0])
        
    return states, progress



def format_dynverse_dataset(*,
    outfile,
    cell_ids,
    tree,
    state, 
    percentage,
    pseudotime,
    mixing_weights,
):

    cell_info = pd.DataFrame(mixing_weights.astype(str), columns=['mix_weight_' + str(i) for i in range(mixing_weights.shape[-1])])
    cell_info['cell_id'] = cell_ids

    branch_progressions = pd.DataFrame(
        {
            'cell_id' : cell_ids, 
            'branch_id' : state.astype(str),
            'percentage' : percentage
        }
    )

    branch_names =  np.unique(state)
    branches_df = pd.DataFrame(
        {'branch_id' : branch_names.astype(str), 
        'directed' : [True]*len(branch_names), 
        'length' : [1]*len(branch_names)}
    )

    edges = np.array(list(nx.dfs_edges(tree, source = -1))[1:]).astype(str)
    branch_network = pd.DataFrame(
        edges, columns = ['from', 'to']
    )

    start_cell = cell_ids[np.argmin(pseudotime)]
    end_cells = branch_progressions[
            branch_progressions.groupby('branch_id')['percentage'].transform(max) == branch_progressions['percentage']
        ].set_index('branch_id').loc[[str(x) for x in tree.nodes() if tree.out_degree(x)==0]].cell_id.values

    trajectory = dyn.wrap_data(cell_ids = cell_ids)\
            .add_branch_trajectory(
                branch_network = branch_network,
                branches = branches_df,
                branch_progressions = branch_progressions,
            ).add_prior_information(
                start_id = [start_cell],
                end_id = end_cells
            ).add_root([start_cell])\
            .add_cell_waypoints()

    trajectory.write_output(outfile)

    write_cell_info(outfile, {
        'mix_weight_' + str(i) : mixing_weights[:, i] for i in range(mixing_weights.shape[1])
    })

    logging.info('Wrote dynverse dataset: ' + outfile)

    
def generate_frankentrajectory(
    state_compositions,
    outfile,
    branch_times = None,
    n_cells = 1000,
    gamma = 0.1,
    max_depth = 2,
    max_width = 2,
    ptime_alpha = 1.,
    ptime_beta = 1.,
    sigmoid_approach = True,
    sigmoid_aggression = 5,
    seed = None,
):

    np.random.seed(seed=seed)
    if branch_times is None:
        branch_times = np.linspace(0, 1., max_depth + 1)
    else:
        branch_times = np.array(branch_times).reshape(-1)
        assert(len(branch_times) == (max_depth + 1) )
    
    state_compositions = np.array(state_compositions)
    assert(len(state_compositions) == max_width**max_depth)

    _, paths, transition_mixing, states, real_time, percentage = \
            generate_branching_process(
                branch_times,
                n_cells=n_cells,
                max_depth=max_depth,
                max_width=max_width,
                gamma = gamma,
                ptime_alpha= ptime_alpha,
                ptime_beta=ptime_beta,
                sigmoid_aggression=sigmoid_aggression,
                sigmoid_approach=sigmoid_approach
            )

    mixing_weights = get_mixing_weights(
            state_compositions = state_compositions, 
            paths = paths, 
            transition_mixing = transition_mixing
    )

    tree = make_graph(states)
    state, percentage = simplify_tree(tree, state_compositions, states, percentage)

    cell_ids = np.array(
        ['Cell_' + str(i) for i in range(len(state))]
    )

    format_dynverse_dataset(
        cell_ids = cell_ids,
        tree = tree,
        state = state, 
        mixing_weights = mixing_weights, 
        percentage = percentage,
        pseudotime = real_time,
        outfile = outfile,
    )


def add_arguments(parser):

    parser.add_argument('--state-composition','-s', type = float, nargs='+', action='append')
    parser.add_argument('--branch-times', '-b', type = float, nargs='+')
    parser.add_argument('--outfile', '-o', type = str, required = True)
    parser.add_argument('--sigmoid-aggression', '-a', type = float, default=5)
    parser.add_argument('--no-sigmoid', action = 'store_const', const = True,
        default = False)
    parser.add_argument('--ptime-beta', '-ptb', type = float, default= 1.)
    parser.add_argument('--ptime-alpha', '-pta', type = float, default=1.)
    parser.add_argument('--max-depth', type = int, default = 2)
    parser.add_argument('--max_width', type = int, default= 2)
    parser.add_argument('--gamma', type = float, default= 0.1)
    parser.add_argument('--n-cells', '-n', type = int, default=1000)
    parser.add_argument('--seed', type = int, default=None)