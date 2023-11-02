#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:15:10 2023

@author: mkw5910

This script calculates the pairwise Fst statistics for all unique pairs of populations in a given tree sequence.
It takes as input a tree sequence file, simulation model name, simulation replicate iteration, and an optional output directory and file.
The script first identifies all unique populations in the tree sequence and then calculates the Fst statistics for each pair of populations.
The results, along with the simulation model name and iteration number, are then written to a tab-separated text file.
"""
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

ap.add_argument("-sn", "--simName", required=True,
   help="Simulation Model Name")

ap.add_argument("-i", "--iteration", required=True,
   help="Simulation replicate iteration")

ap.add_argument("-o", "--out", required=False,
   help="Output Directory and File")

args = vars(ap.parse_args())

# Define variables
iteration = args['iteration']
out = args['out']
sim_model = args['simName']


## **** ### ***** #### ## **** ### ***** #### ## **** ### ***** ####
# read in the tree-sequence
ts_branch = tskit.load(args["tree_sequence"])

# restrict to only the nodes that have simulated samples
sample_nodes = ts_branch.tables.nodes[ts_branch.tables.nodes.flags==1]

# Get the unique list of population IDs (0,1,2,..)
unique_populations = np.unique(sample_nodes.population)

# Generate an empty dictionary with the population labels as the entry and pop order [0,1,2,..] as the key
pop_node_Dict = {}
for i in range(len(unique_populations)):
        idx = unique_populations[i]
        #print(i)
        #print(ts_branch.tables.populations[i])
        #print(ts_branch.tables.populations[i].metadata['name'])
        pop_node_Dict[idx] = [ts_branch.tables.populations[idx].metadata['name']]
        
# Add the ind nodes to the dict for each population key
for i in range(len(unique_populations)):
        #print(unique_populations[i])
        idx = unique_populations[i]
        pop_inds = ts_branch.samples(population=idx)
        pop_node_Dict[idx].append(pop_inds)

# Generate all unique 2-way combinations of populations 
pop_combinations = list(itertools.combinations(unique_populations, 2))

# Create an empty list to store the dictionaries [pop1, pop2, fst, sim iteration & sim model]
data = []
for i in range(len(pop_combinations)):
    #print(pop_combinations[i])
    pop1 = pop_combinations[i][0]
    pop2 = pop_combinations[i][1]
    fst = ts_branch.Fst(sample_sets=[pop_node_Dict[pop1][1].tolist(), pop_node_Dict[pop2][1].tolist()], windows = [0.0, ts_branch.sequence_length], mode='branch', span_normalise=True)
    
    # Create a dictionary for the current row
    row = {
        'pop1': pop_node_Dict[pop1][0],
        'pop2': pop_node_Dict[pop2][0],
        'fst': fst[0],
        'iteration': iteration,
        'sim_model': sim_model
    }
    
    # Append the dictionary to the list
    data.append(row)

# Create a dataframe from the list of dictionaries
df = pd.DataFrame(data)

# write dataframe to file as tab-separated txt file
df.to_csv(f"{out}", sep='\t', index=False)
## **** ### ***** #### ## **** ### ***** #### ## **** ### ***** ####
## **** ### ***** #### ## **** ### ***** #### ## **** ### ***** ####








