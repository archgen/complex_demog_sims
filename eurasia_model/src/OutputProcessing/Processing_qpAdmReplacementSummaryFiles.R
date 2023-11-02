#' Reading in the complex simulation qpAdm Replacement aDNA and branch
#' files, merging them and generating a summary analysis file

#' Load libraries.
library(data.table)

#' set working directory to output top directory
setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/out/complexDemography_model_A_HapMapII_GRCh37__Chrs21_22/")

#' qpAdm Replacement File
branch_qpAdmReplace = fread(file = "./admixtools2/ts_branch_SummaryOutput_qpAdmReplacement/ts_branch_qpAdmReplacement_SummaryTable__complexDemography_model_A_HapMapII_GRCh37__Chrs21_22.txt")
aDNA_qpAdmReplace = fread(file = "./admixtools2/ts_mu_eig_aDNA_SummaryOutput_qpAdmReplacement/ts_mu_eig_aDNA_qpAdmReplacement_SummaryTable__complexDemography_model_A_HapMapII_GRCh37__Chrs21_22.txt")
