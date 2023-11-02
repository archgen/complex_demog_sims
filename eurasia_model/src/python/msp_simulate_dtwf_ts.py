#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 11:46:27 2022

@author: mkw5910
"""
import sys
import msprime
print('msprime', msprime.__version__)
import pandas as pd
from tabulate import tabulate
import demes
import argparse
import random
import math

# SYSTEM ARGUMENTS

# E.G  https://www.datacamp.com/tutorial/argument-parsing-in-python
ap = argparse.ArgumentParser()
ap.add_argument("--simName", required=True,
   help="Simulation Model Name")

ap.add_argument("--demes", required=True,
   help="Path to Demes yaml file")

ap.add_argument("--samples", required=True,
   help="Path to sample sheet .txt file")

ap.add_argument("--WFgens", required=True,
   help="Number of generations to simulate under Wright Fisher")

ap.add_argument("--iteration", required=True,
   help="Simulation replicate iteration")

ap.add_argument("--out", required=True,
   help="Output Directory and File")

args = vars(ap.parse_args())


# SIMULATION PARAMETERS
#######################
sim_model = args['simName']
deme_yaml = demes.load(args['demes'])
sample_sheet = pd.read_csv(args['samples'])
DTWFgens = format(int(float(args['WFgens'])))
iteration = format(int(float(args['iteration'])))
out = args['out']
random.seed()
seed = random.randint(1, 2**32)

# Samples to simulate
samp_Ns = list(sample_sheet.Nsamples)
pop_list = list(sample_sheet.Pop)
samp_gens = list(sample_sheet.Generations)
sList = list(zip(samp_Ns,pop_list,samp_gens))
samples_list = []
for i in sList: samples_list.append(msprime.SampleSet(num_samples=i[0], population=i[1], time=i[2]))
print(tabulate(sample_sheet, headers=["Number of Samples", "Demes Population", "Generation Sampled", "qpAdm Rotation Analysis", "Output Format", "qpAdm Fixed Analysis"]))

# Genome Recombination Map
chr1_L = 248956422
chr2_L = 242193529
chr3_L = 198295559
chr4_L = 190214555
chr5_L = 181538259
chr6_L = 170805979
chr7_L = 159345973
chr8_L = 145138636
chr9_L = 138394717
chr10_L = 133797422
chr11_L = 135086622
chr12_L = 133275309
chr13_L = 114364328
chr14_L = 107043718
chr15_L = 101991189
chr16_L = 90338345
chr17_L = 83257441
chr18_L = 80373285
chr19_L = 58617616
chr20_L = 64444167
chr21_L = 46709983
chr22_L = 50818468


chr1_r = 1.15235e-08
chr2_r = 1.10429e-08
chr3_r = 1.12585e-08
chr4_r = 1.1482e-08
chr5_r = 1.12443e-08
chr6_r = 1.12659e-08
chr7_r = 1.17713e-08
chr8_r = 1.16049e-08
chr9_r = 1.21987e-08
chr10_r = 1.33337e-08
chr11_r = 1.17213e-08
chr12_r = 1.30981e-08
chr13_r = 1.3061e-08
chr14_r = 1.36298e-08
chr15_r = 1.73876e-08
chr16_r = 1.48315e-08
chr17_r = 1.55383e-08
chr18_r = 1.46455e-08
chr19_r = 1.83848e-08
chr20_r = 1.67886e-08
chr21_r = 1.72443e-08
chr22_r = 2.10572e-08


#' formaing WH Map
chrom_lengths = [chr1_L, chr2_L, chr3_L, chr4_L, chr5_L,
chr6_L, chr7_L, chr8_L, chr9_L, chr10_L, chr11_L, chr12_L,
chr13_L, chr14_L, chr15_L, chr16_L, chr17_L, chr18_L, chr19_L,
chr20_L, chr21_L, chr22_L]
chrom_positions = [0]
for i in range(len(chrom_lengths)):
    chrom_positions.append(chrom_positions[i] + chrom_lengths[i])

map_positions = [chrom_positions[0]]
for position in chrom_positions[1:]:
    map_positions.extend([position, position + 1])
map_positions = map_positions[0:-1]


r_break = math.log(2)
chr_rates = [chr1_r, chr2_r, chr3_r, chr4_r, chr5_r,
chr6_r, chr7_r, chr8_r, chr9_r, chr10_r, chr11_r, chr12_r,
chr13_r, chr14_r, chr15_r, chr16_r, chr17_r, chr18_r, chr19_r,
chr20_r, chr21_r, chr22_r]
rates = []
for rate in chr_rates:
    rates.extend([rate, r_break])
rates = rates[0:-1]


rate_map = msprime.RateMap(position=map_positions, rate=rates)

# Demography
demography = msprime.Demography.from_demes(deme_yaml)
#print(demography.debug())

print("")
print("RECOMBINATION MAP")
print(rate_map)



print(f" Simulation of python msp simulation model {sim_model}. DTWF gennerations: {DTWFgens} replicate: {iteration} seed: {seed}")
# MSPRIME DTWF SIMULATION
ts = msprime.sim_ancestry(
    samples=samples_list,
    demography=demography,
    recombination_rate=rate_map,
    model = [
        msprime.DiscreteTimeWrightFisher(duration=int(DTWFgens)),
        msprime.StandardCoalescent(),],
    random_seed=seed,
    record_provenance=False)
print(ts)
ts.dump(f"{out}")









