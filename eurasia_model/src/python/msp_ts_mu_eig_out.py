#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Oct 24 13:38:39 2022

@author: mkw5910
"""
import os
import msprime
print('msprime', msprime.__version__)
import argparse
import tskit
import pandas as pd

# SYSTEM ARGUMENTS
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-sn", "--simName", required=True,
   help="Simulation Model Name")
ap.add_argument("-ts", "--tree_sequence", required=True,
   help="Path to tree sequence")
ap.add_argument("-i", "--iteration", required=True,
   help="Simulation replicate iteration")
ap.add_argument("-ss", "--samples", required=True,
   help="Path to sample sheet .txt file")
ap.add_argument("-o", "--out", required=True,
   help="Output Directory and File")
args = vars(ap.parse_args())

# INPUT PARAMETERS
#######################
sample_sheet = pd.read_csv(args['samples'])
ts_mu = tskit.load(args["tree_sequence"])
sim_model = args['simName']


# Generating sample information
# Creating sample table with # sampels, population, generation and their output status 
samp_Ns = list(sample_sheet.Nsamples)
pop_list = list(sample_sheet.Pop)
sampTime = list(sample_sheet.Generations)
output = list(sample_sheet.Output)
sList = list(zip(samp_Ns,pop_list,sampTime,output))
samples_list = []
for i in sList: samples_list.append(msprime.SampleSet(num_samples=i[0], population=i[1], time=i[2]))

# Creating a list of popualtions for all the simulated sampled populations that are marked as eigenstrat out
popualtion_eig = list()
for row in range(len(samples_list)):
    if sList[row][3] == "eig":
        popualtion_eig.append(f'{samples_list[row].population}')
# Retrieving the index position in the genotyoe files for all samples from populations that have eigenstrat set as output
samp_eig_index = list()
samp_eig_inds = list()
pop_eig = list()
for row in range(ts_mu.sample_size):
    if f'{ts_mu.tables.populations[ts_mu.tables.nodes[row].population].metadata["name"]}' in popualtion_eig:
        samp_eig_index.append(row)
        samp_eig_inds.append(f'ind_{row}_{ts_mu.tables.populations[ts_mu.tables.nodes[row].population].metadata["name"]}__gen_{round(ts_mu.tables.nodes[row].time)}')
        pop_eig.append(f'{ts_mu.tables.populations[ts_mu.tables.nodes[row].population].metadata["name"]}')


# Generating the .geno and .snp files
# Set output DIR
os.chdir(args["out"])
chrom=1
write_eig_geno = list()
write_eig_snp = list()
for variant in ts_mu.variants(): 
    geno_lst = str(variant.genotypes[samp_eig_index].tolist())
    geno_str = str(geno_lst)[1:-1].replace(",","").replace(" ", "")
    Eig_geno = geno_str.replace("0","2").replace("1","0")
    write_eig_geno.append(Eig_geno)
    write_eig_snp.append(f'ss_{variant.index}\t{chrom}\t{variant.position/float(ts_mu.sequence_length)}\t{variant.position}\t2\t0')
 # Write output files    
with open('ts_mu_eig_%s_rep_{}.geno'.format(args['iteration']) %(args['simName']), mode='wt', encoding='utf-8') as geno:
    geno.write('\n'.join(write_eig_geno)) # write .geno file not in the loop as large file
with open('ts_mu_eig_%s_rep_{}.snp'.format(args['iteration']) %(args['simName']), mode='wt', encoding='utf-8') as snp:
    snp.write('\n'.join(write_eig_snp)) # write .snp file not in the loop as large file


with open('ts_mu_eig_%s_rep_{}.ind'.format(args['iteration']) %(args['simName']), 'w') as indFile:
    for s, p in zip(samp_eig_inds, pop_eig):
        indFile.write("%s\tU\t%s \n" %(s,p))
    indFile.close()
