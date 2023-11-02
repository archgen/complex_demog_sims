# Code for Testing Times: Admixture Inference Challenges 

This repository contains the code for reproducing simulations from Williams et al. (2023), Testing Times: Challenges in Disentangling Admixture Histories in Recent and Complex Demographies.


## Overview
Below are the two simulation project sub-directories (simple_models, and eurasia_model) that each contain Snakemake pipelines to execute them.

- simple_models: contains the Snakemake pipeline to simulate the simple demographic model simulations (Model 1-4) that explore the specified topological and admixture parameter space.
- eurasian_model: contains the Snakemake pipeline to simulate from the Eurasian complex demographic model.  

The Snakemake pipelines were executed from a (snake_base) conda environment using anaconda Command line client (version 1.7.2), Snakemake v.3.13.3, and Python 3.6.13 :: Anaconda, Inc. When running the pipeline for the first time, the Snakemake software will initialize a conda environment for the simulation pipeline. The file, *sims-env.yaml* contains all the software necessary to execute the simulations and Snakemake will read from this file to set-up the conda env (note conda can be finicky when transferring across operating systems, the pipeline was setup on a Red Hat Enterprise Linux 7 env.). 


## simple_models
- src: contains the python, and r code. 
- config: contains the simulation configuration files for each demographic model, and a template cluster configuration file (slurm).  

The simple_models root directory contains the executable Snakemake pipeline smk script *snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk* which you execute through

`snakemake -s snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk --configfile config/config_msp_dtwf_modelX.yml --use-conda`

You can override the configuration file through the command line argument `--config yourparam=1.5.` As such to perform a small chromosome test run the following:

`snakemake -s snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk --configfile config/config_msp_dtwf_model[1-4].yml --use-conda --config chromosome_length=1e6 recombination_rate=1e-8 window_length_Fstats=1e5`


## eurasia_model
- src: contains the bash, python, and r code. 
- config: contains the simulation parameters configuration code.  
- input: contains the demes demography, sample sheet.
- dev_RLibs: contains the development version of data.table for fast writing of Eigenstrat files.
