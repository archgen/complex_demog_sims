#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:04:36 2022

@author: mkw5910
"""
import msprime
print('msprime', msprime.__version__)
import argparse
import tskit
import random
from tabulate import tabulate


# SYSTEM ARGUMENTS
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-ts", "--tree_sequence", required=True,
   help="Path to tree sequence")

ap.add_argument("-mu", "--mutation_rate", required=True,
   help="Mutation rate")

ap.add_argument("--WFgens", required=True,
   help="Number of generations to simulate under Wright Fisher")

ap.add_argument("--length", required=True,
   help="Indicator of genome length (chr22 or WG)")

ap.add_argument("-i", "--iteration", required=True,
   help="Simulation replicate iteration")

ap.add_argument("-o", "--out", required=True,
   help="Output Directory and File")
args = vars(ap.parse_args())

# INPUT PARAMETERS
#######################
out = args["out"]
# Read in msprime tree-sequence 
ts = tskit.load(args["tree_sequence"])
random.seed()
seed = random.randint(1, 2**32)
DTWFgens = format(int(float(args['WFgens'])))
iteration = format(int(float(args['iteration'])))
length = format(args['length'])

ts_mu = msprime.sim_mutations(ts, rate=args['mutation_rate'], random_seed=seed)
print(f'Simulated {ts_mu.num_mutations} mutations with a mutation rate of {args["mutation_rate"]}')
ts_mu.dump(f"{out}")

# Output simulation paramters table
table = [["DTWFgens", "rep", "seed", "ts_length", "N_ts_mu_trees", "N_ts_mu_nodes", "N_ts_mu_edges", "N_ts_mu_sites", "N_ts_mu_mutations", "mutation_rate", "genome_length"],
                 [DTWFgens, iteration, seed, ts.sequence_length, ts_mu.num_trees, ts_mu.num_nodes, ts_mu.num_edges, ts_mu.num_sites, ts_mu.num_mutations, args['mutation_rate'], length]]
print(tabulate(table))
#print("Source sampling at generation time = ", (T_admix / generation_time) + generation_time)
with open(f"{out}_ts_mu_summary.txt", 'w') as f:
    f.write(tabulate(table))
