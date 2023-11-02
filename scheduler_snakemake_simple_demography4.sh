#!/bin/bash
#SBATCH --partition=open 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --time=48:00:00 
##SBATCH --time=01:00:00
#SBATCH --mem=128GB 
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mkw5910@psu.edu


cd /storage/group/cdh5313/default/mkw5910/sims/github

mkdir -p logs_slurm
snakemake -s snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk --configfile config_model2.yml --use-conda --cluster-config cluster.yaml --cluster 'sbatch -t {cluster.time} --mem {cluster.mem} -c {cluster.cpus} --cores 1 -N 1 -o {cluster.output} -e {cluster.output} --mail-type {cluster.email_type} --mail-user {cluster.email}' -j 99 --latency-wait 604800 --restart-times 0 --keep-going
