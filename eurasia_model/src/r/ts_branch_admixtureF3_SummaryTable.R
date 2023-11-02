#' Generate admixture F3 Summary Table
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(tidyverse)
library(optparse)
library(data.table)

# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("-d", "--in_dir"),
              type = "character",
              default = NULL,
              help = "In Dir for qpAdm RDS files",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "output FILENAME & PATH",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
arguments = parse_args(opt_parser, positional_arguments=TRUE);
opt <- arguments$options

if (is.null(opt$in_dir)){
  print_help(opt_parser)
  stop("In DIR for admixture F3 RDS files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}

# message(" ****** DEBUGGING INPUT PARAMS ****** ")
# # opt_parser = OptionParser(option_list = option_list);
# # opt = parse_args(opt_parser);
# message("{params.workDIR} = ", opt$in_dir)
# message("{params.datatype} = ", opt$data_type)
# message("{params.sim_model} = ", opt$model_name)
# message("{params.qpAdmTable} = ", opt$out)
# message(" ****** END  ****** ")


# Parameters
#' Creating a list of the rds files for all the replicates for each target
admixF3_reps_SummaryTable = data.table()
in_rds_list = list.files(opt$in_dir)
print(in_rds_list)
for(i in 1:length(in_rds_list)){
  admixF3_rep = data.table(readRDS(file = paste0(opt$in_dir, in_rds_list[i])))
  #admixF3_rep = data.table(readRDS(file = paste0(in_dir, in_rds_list[i]))) ##! TESTING
  admixF3_DT = admixF3_rep[!duplicated(admixF3_rep[, c("est", "se", "z", "p")]), ]
  admixF3_reps_SummaryTable = rbind(admixF3_reps_SummaryTable, admixF3_DT)
}

# Save data table as RDS
saveRDS(admixF3_reps_SummaryTable, file=opt$out)

#' ####################################################################################################################################################################################################
#' TESTING
# setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/project__complex_msp_simulation/")
#in_dir = "out/complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource/admixtools2/ts_branch_admixtureF3/"
#in_rds_list = list.files(in_dir)
# in_prefix = "out/complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource/admixtools2/admixtureF3_aDNA/admixtureF3_sim_aDNA__complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource_rep_"
# admixF3_rep = readRDS(file = paste0(inDIR, in_rds_list[1]))




