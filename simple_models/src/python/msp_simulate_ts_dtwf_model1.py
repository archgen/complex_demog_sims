#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 13:09:27 2022

Simple seven population demography. 

@author: mkw5910

This script simulates a simple seven population demography using msprime. It takes as input various parameters such as the simulation model name, 
simulation replicate iteration, generations to simulate under the DTWF model, generation time for the deepest population split (T0), simulation 
generation time, simulation diploid effective population size, number of diploid samples to simulate, and an output directory and file. The script 
first randomly samples simulation parameters until certain conditions are met. Then it sets up a demography with seven populations and simulates 
their evolution according to the specified parameters. The results are then written to a file.
"""

# IMPORT PACKAGES
import msprime
from tabulate import tabulate
import argparse
import random
import numpy as np
import math

print('msprime', msprime.__version__)
print('argparse', argparse.__version__)
print('np', np.__version__)

# SYSTEM ARGUMENTS
#######################
ap = argparse.ArgumentParser()
ap.add_argument("--simName", required=True,
   help="Simulation Model Name")

ap.add_argument("--length", required=False, default=None,
   help="Single-genome simulation sequence length")

ap.add_argument("--recom", required=False, default=None,
   help="Single-genome simulation sequence recombination rate")

ap.add_argument("--iteration", required=True,
   help="Simulation replicate iteration")

ap.add_argument("--DTWFgens", required=True,
   help="Generations to simulate under DTWF model")

ap.add_argument("--MRCAgen", required=True,
   help="Generation time for the deepest population split (T0)")

ap.add_argument("--genTime", required=True,
   help="Simulation generation time")

ap.add_argument("--Ne", required=True,
   help="Simulation diploid effecitve population size")

ap.add_argument("--Nsamples", required=True,
   help="Number of diploid samples to simulate")

ap.add_argument("--out", required=True,
   help="Output Directory and File")

args = vars(ap.parse_args())

# Parameters Defining
sim_model = args['simName']
iteration = int(float(args['iteration']))
DTWFgens = int(float(args['DTWFgens']))
MRCAgen = int(float(args['MRCAgen']))
Nsamples = int(float(args['Nsamples']))
out = args['out']
random.seed()
seed = random.randint(1, 2**32)
Ne_base = int(float(args['Ne']))
generation_time = int(float(args['genTime']))


# Random Simulation Parameter Sampling 
while True:
    T_0 = MRCAgen
    T_1 = int(np.random.uniform(generation_time * 4, T_0, 1))
    T_2 = int(np.random.uniform(generation_time * 3, T_1, 1))
    T_3 = int(np.random.uniform(generation_time * 2, T_2, 1))
    T_4 = int(np.random.uniform(generation_time * 2, T_2, 1))
    T_admix = int(np.random.uniform(generation_time, min(T_3, T_4), 1))
    alpha = np.random.uniform(0, 1, 1)[0]
    
    if T_1 > T_2 > T_3 > T_admix and T_1 > T_2 > T_4 > T_admix:
        break


# Demography Parameters
# Output simulation paramters table
sim_params_table = [["T_0_gen", "T_1_gen", "T_2_gen", "T_3_gen", "T_4_gen", "T_admix_gen", "alpha_S1", "DTWFgens", "generation_time", "rep", "seed"],
                 [T_0 / generation_time, T_1 / generation_time, T_2 / generation_time, T_3 / generation_time, T_4 / generation_time, T_admix / generation_time, alpha, DTWFgens, generation_time, iteration, seed]]
print(tabulate(sim_params_table))
    
# Demography
demography = msprime.Demography()
demography.add_population(name="R4", initial_size=Ne_base, initially_active=True)
demography.add_population(name="R3", initial_size=Ne_base, initially_active=True) 
demography.add_population(name="R2", initial_size=Ne_base, initially_active=True)
demography.add_population(name="R1", initial_size=Ne_base, initially_active=True)
demography.add_population(name="S2", initial_size=Ne_base, initially_active=True)
demography.add_population(name="S1", initial_size=Ne_base, initially_active=True)
demography.add_population(name="T", initial_size=Ne_base, initially_active=True)
demography.add_population(name="R2S2", initial_size=Ne_base)
demography.add_population(name="R1S1", initial_size=Ne_base)
demography.add_population(name="iR12iS12", initial_size=Ne_base)
demography.add_population(name="iANC", initial_size=Ne_base)
demography.add_population(name="ANC", initial_size=Ne_base, initially_active=True)


# Demography Simulation Set-up
# ADMIXTURE: T as mixture of S1 and S2. 
demography.add_admixture(time=T_admix / generation_time, derived="T", ancestral=["S1", "S2"], proportions=[alpha, 1-alpha])

# SPLIT: R2 and S2 from ancestral S2-R2
demography.add_population_split(time=T_4 / generation_time, derived=["R2", "S2"], ancestral="R2S2") 

# SPLIT: R1 and S1 from ancestral S1-R1
demography.add_population_split(time=T_3 / generation_time, derived=["R1", "S1"], ancestral="R1S1") 

# SPLIT: Ancestral R1-S1 and R2-S2 from internal ancestral S12-R12
demography.add_population_split(time=T_2 / generation_time, derived=["R1S1", "R2S2"], ancestral="iR12iS12") 

# SPLIT: R3 and internal S12-R12 from ancestral internal
demography.add_population_split(time=T_1 / generation_time, derived=["R3", "iR12iS12"], ancestral="iANC")

# SPLIT: internal ancestral and R4 from ancestral 
demography.add_population_split(time=T_0 / generation_time, derived=["R4", "iANC"], ancestral="ANC")

# Sort demography events as T_4 may be larger than T_3
demography.sort_events()

# Tree-sequence Simulation for single genome/chromosome 
if args['recom'] and args['length'] != "None":
        seq_L = int(float(args['length']))
        seq_r = float(args['recom'])
        print("Sequence length is ", seq_L)
        print("Recombination rate is ", seq_r)

        # MSPRIME SIMULATIONS
        ts = msprime.sim_ancestry(
            samples = [
                msprime.SampleSet(num_samples = Nsamples, population="ANC", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="R4", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="R3", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="R2", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="R1", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="S2", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="S1", time=0),
                msprime.SampleSet(num_samples = Nsamples, population="T", time=0)
            ],
            #samples=samples_list,
            demography=demography,
            recombination_rate=seq_r,
            sequence_length=seq_L,
            model = [
                msprime.DiscreteTimeWrightFisher(duration=DTWFgens),
                msprime.StandardCoalescent(),],
            random_seed=seed)
        print(ts)
        # Output tree-sequence
        ts.dump(f"{out}")
else:
    chr1 = 248956422
    chr2 = 242193529
    r_chrom = [1.15235e-08, 1.10429e-08]
    r_break = math.log(2)
    chrom_positions = [0, chr1, (chr1+chr2)]
    map_positions = [
        chrom_positions[0],
        chrom_positions[1],
        chrom_positions[1] + 1,
        chrom_positions[2]]
    rates = [r_chrom[0], r_break, r_chrom[1]]
    rate_map = msprime.RateMap(position=map_positions, rate=rates)
    print("Rate map popsim-HomSap Chromosomes 1-2 ", rate_map)

    # MSPRIME SIMULATIONS
    ts = msprime.sim_ancestry(
        samples = [
            msprime.SampleSet(num_samples = Nsamples, population="ANC", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="R4", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="R3", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="R2", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="R1", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="S2", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="S1", time=0),
            msprime.SampleSet(num_samples = Nsamples, population="T", time=0)
        ],
        demography=demography,
        recombination_rate=rate_map,
        model = [
            msprime.DiscreteTimeWrightFisher(duration=DTWFgens),
            msprime.StandardCoalescent(),],
        random_seed=seed)
    print(ts)
    # Output tree-sequence
    ts.dump(f"{out}")

# Output simulation paramters table
table = [["T_0_gen", "T_1_gen", "T_2_gen", "T_3_gen", "T_4_gen", "T_admix_gen", "alpha_S1",  "DTWFgens", "generation_time", "rep", "seed", "N_ts_trees", "N_ts_nodes", "N_ts_edges"],
                 [T_0 / generation_time, T_1 / generation_time, T_2 / generation_time, T_3 / generation_time, T_4 / generation_time, T_admix / generation_time, alpha, DTWFgens, generation_time, iteration, seed, ts.num_trees, ts.num_nodes, ts.num_edges]]
print(tabulate(table))
#print("Source sampling at generation time = ", (T_admix / generation_time) + generation_time)
with open(f"{out}_simulation_params.txt", 'w') as f:
    f.write(tabulate(table))
