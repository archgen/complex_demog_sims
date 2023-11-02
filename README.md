# Code for Testing Times: Admixture Inference Challenges 

This repository contains the code for reproducing simulations from Williams et al. (2023), Testing Times: Challenges in Disentangling Admixture Histories in Recent and Complex Demographies.


## Overview
Below are the two simulation project sub-directories (simple_models, and eurasia_model) that each contain Snakemake pipelines to execute them.

- simple_models: contains the Snakemake pipeline to simulate the simple demographic model simulations (Model 1-4) that explore the specified topological and admixture parameter space.
- eurasian_model: contains the Snakemake pipeline to simulate from the Eurasian complex demographic model.  

The Snakemake pipelines were executed from a (snake_base) conda environment using anaconda Command line client (version 1.7.2), Snakemake v.3.13.3, and Python 3.6.13 :: Anaconda, Inc. When running the pipeline for the first time, the Snakemake software will initialize a conda environment for the simulation pipeline. The file, *msprime-env.yaml*, located in each sub-directory, contains all the software necessary to execute the simulations and Snakemake will read from this file to set-up the conda env (note conda can be finicky when transferring across operating systems, the pipeline was setup on a Red Hat Enterprise Linux 7 env.). The only other software required is the admixools2 R package (https://uqrmaie1.github.io/admixtools/articles/admixtools.html).


## simple_models
- src: contains the python, and r code. 
- config: contains the simulation configuration files for each demographic model.  

The simple_models root directory contains the executable Snakemake pipeline smk script *snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk* which you execute through

`snakemake -s snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk --configfile config/config_msp_dtwf_modelX.yml --use-conda`

You can override the configuration file through the command line argument `--config yourparam=1.5.` As such to perform a small chromosome test run the following from the terminal:

`snakemake -s snakemake_msp_dtwf__qpAdm_f3_SimpleDemog.smk --configfile config/config_msp_dtwf_model[1-4].yml --use-conda --config chromosome_length=1e6 recombination_rate=1e-8 window_length_Fstats=1e5`


## eurasia_model
- src: contains the bash, python, and r code. 
- config: contains the simulation parameters configuration code.  
- input: contains the demes demography, sample sheet, and where to store the eig aDNA datasets (see below).
- dev_RLibs: contains data.table v.1.14.9 dev version to access the fwrite function, and gtools Version: 3.9.4. 

The eurasia_model root directory contains the executable Snakemake pipeline smk script *snakemake_msp_dtwf_ts_mu_eigaDNA__eurasia_qpAdm_f3.smk* which you execute through
`snakemake -s snakemake_msp_dtwf_ts_mu_eigaDNA__eurasia_qpAdm_f3.smk --configfile config/config_msp_dtwf_eurasia.yml --use-conda 

#### Empirical data misssingness 
In order to perform the data missingness aDNA step please download the three aDNA datasets to the input directory which can be accessed at the following link: [eig_aDNA_data_dir](https://drive.google.com/drive/folders/1Uv-2NSK7e-EtO960sKGkEHB_xV5k9bOl?usp=sharing)

