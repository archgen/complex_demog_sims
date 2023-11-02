#' Generate qpAdm Summary Table for the qpAdm Replacement analyses
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(tidyverse)
library(optparse)
#library(gtools, lib.loc="/storage/group/cdh5313/default/mkw5910/sims/complex_Demography/dev_RLibs")
library(gtools)
library(data.table)

##### TESTING
#setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/out/complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource/")
#model = "snake_AncientEurasia_9K19"
#inDIR = "./admixtools2/ts_branch_qpAdmReplacement_Table/"
#results_files_list = mixedsort(sort(list.files(inDIR, full.names = TRUE)))

# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("-d", "--in_dir"),
              type = "character",
              default = NULL,
              help = "In Dir for qpAdm RDS files",
              metavar = "character"),
  make_option(c("--data_type"),
              type = "character",
              default = NULL,
              help = "Data type used to generate qpAdm analysis",
              metavar = "character"),
  make_option(c("-m", "--model_name"),
              type = "character",
              default = NULL,
              help = "Name of the simulation Model",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "output FILENAME & PATH",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$in_dir)){
  print_help(opt_parser)
  stop("In DIR for qpAdm RDS files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}


# Parameters
dataType= opt$data_type
modelName = opt$model_name

# Read in summary tables 
results_files_list = mixedsort(sort(list.files(opt$in_dir, full.names = TRUE)))
qpadm_RDS <- lapply(results_files_list, readRDS)

# Generate Summary Table
qpAdm_fixed_summary_table = data.table()
for(i in 1:length(qpadm_RDS)){
  qpAdmDT = data.frame(qpadm_RDS[i])
  qpAdmDT = data.table(qpAdmDT)
  qpAdmDT[, dataType:=dataType]
  qpAdmDT[, simModelName:=modelName]
  qpAdmDT[, summary_table_rep:=i]
  qpAdmDT[, plausible := (weight.p1 >= 0 & weight.p1 <= 1 & weight.p2 >=0 & weight.p2 <=1 & p.value >= 0.01)]
  qpAdm_fixed_summary_table = rbind(qpAdm_fixed_summary_table, qpAdmDT)
}

# Write output Table
write.table(qpAdm_fixed_summary_table, file=opt$out, row.names = FALSE)



