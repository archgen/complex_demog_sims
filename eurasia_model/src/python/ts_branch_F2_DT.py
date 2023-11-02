#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:15:10 2023

@author: mkw5910
"""

import sys

# Add the directory to sys.path
#sys.path.insert(0, "/storage/home/mkw5910/.conda/envs/msprime-env/lib/python3.9/site-packages/")

# IMPORT PACKAGES
import msprime
import pandas as pd
import argparse
import numpy as np
import itertools
import tskit

print('msprime', msprime.__version__)
print('pd', pd.__version__)
print('argparse', argparse.__version__)
print('np', np.__version__)
print('tskit', tskit.__version__)


# SYSTEM ARGUMENTS
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-ts", "--tree_sequence", required=True,
   help="Path to tree sequence")

ap.add_argument("-wS", "--window_size", required=True,
   help="Size of windows to compute f2")

ap.add_argument("-o", "--out", required=False,
   help="Output Directory and File")

args = vars(ap.parse_args())


# Define variables
out = args['out']
window_size = int(float(args['window_size']))

# # Testing variables
#out = "TEST_F2_array.txt"
# ts_branch = tskit.load("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/project__complex_msp_simulation/out/HomSap__AncientEurasia_9K19_HapMapII_GRCh37__Chrs21_22/ts_mu/ts_mu_HomSap__AncientEurasia_9K19_HapMapII_GRCh37__Chrs21_22_rep_1.ts")
#window_size = 5000000
#window_size = int(ts_branch.sequence_length / 4)




## **** ### ***** #### ## **** ### ***** #### ## **** ### ***** ####
# read in the tree-sequence
ts_branch = tskit.load(args["tree_sequence"])

# restrict to only the nodes that have simulated samples
sample_nodes = ts_branch.tables.nodes[ts_branch.tables.nodes.flags==1]

# Get the unique list of population IDs (0,1,2,..)
unique_populations = np.unique(sample_nodes.population)
print("unique_populations  ", unique_populations)
# Generate an empty dictionary with the population labels as the entry and pop order [0,1,2,..] as the key
PopNames = []
pop_node_Dict = {}
for i in range(len(unique_populations)):
        #print(unique_populations[i])
        idx = unique_populations[i]
        #print(ts_branch.tables.populations[i])
        #print(ts_branch.tables.populations[i].metadata['name'])
        pop_node_Dict[idx] = [ts_branch.tables.populations[idx].metadata['name']]
        PopNames.append(ts_branch.tables.populations[idx].metadata['name'])

# Add the ind nodes to the dict for each population key
for i in range(len(unique_populations)):
        #print(unique_populations[i])
        idx = unique_populations[i]
        pop_inds = ts_branch.samples(population=idx)
        pop_node_Dict[idx].append(pop_inds)

print("pop_node_Dict     ", pop_node_Dict)

# Generate all unique 2-way combinations of populations 
pop_combinations = list(itertools.combinations(unique_populations, 2))


# Create the blocks list to compute the F2s
num_windows = int(ts_branch.sequence_length /  window_size)

# Create an empty list to store the dictionaries [pop1, pop2, fst, sim iteration & sim model]
f2_table_window = pd.DataFrame()
for i in range(len(pop_combinations)):
    idx = pop_combinations[i]
    pop1 = pop_combinations[i][0]
    pop2 = pop_combinations[i][1]
    print("pop_node_Dict populations  ", pop_node_Dict[pop1], "  ::  " , pop_node_Dict[pop2])
    f2 = ts_branch.f2(sample_sets=[pop_node_Dict[pop1][1].tolist(), pop_node_Dict[pop2][1].tolist()], windows=np.linspace(0, int(ts_branch.sequence_length), num_windows + 1), span_normalise=True, mode="branch")
    name_pop1 = pop_node_Dict[pop1][0]
    name_pop2 = pop_node_Dict[pop2][0]
    pop_pair = name_pop1 + "_AND_" + name_pop2
    # Create table of f2 values and the pops as column name
    f2_run = pd.DataFrame(data=f2, columns=[pop_pair])

    # Concatenate the pairwise F2s by columns
    f2_table_window = pd.concat([f2_table_window, f2_run], axis=1)

# write table to file as tab-separated txt file
f2_table_window.to_csv(f"{out}", sep='\t', index=False)







# Create the blocks list to compute the F2s
#blocks_list = list(range(0, int(ts_branch.sequence_length), window_size))
#blocks_list.append(ts_branch.sequence_length)

# Create an empty list to store the dictionaries [pop1, pop2, fst, sim iteration & sim model]
#f2_table = pd.DataFrame()
#for i in range(len(pop_combinations)):
#    idx = pop_combinations[i]
#    pop1 = pop_combinations[i][0]
#    pop2 = pop_combinations[i][1]
#    print("pop_node_Dict populations  ", pop_node_Dict[pop1], "  ::  " , pop_node_Dict[pop2])
#    f2 = ts_branch.f2(sample_sets=[pop_node_Dict[pop1][1].tolist(), pop_node_Dict[pop2][1].tolist()], windows = blocks_list, mode='branch', span_normalise=True)
#    name_pop1 = pop_node_Dict[pop1][0]
#    name_pop2 = pop_node_Dict[pop2][0]
#    pop_pair = name_pop1 + "_AND_" + name_pop2
#    # Create table of f2 values and the pops as column name
#    f2_run = pd.DataFrame(data=f2, columns=[pop_pair])
#
#    # Concatenate the pairwise F2s by columns
#    f2_table = pd.concat([f2_table, f2_run], axis=1)
#
## write table to file as tab-separated txt file
#f2_table.to_csv(f"{out}", sep='\t', index=False)
