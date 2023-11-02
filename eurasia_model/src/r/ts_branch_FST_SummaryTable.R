#' Generate branch average FST Summary Table
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
#.libPaths("/storage/group/cdh5313/default/mkw5910/sims/complex_Demography/.snakemake/conda/5af5b266/lib/R/library/")
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
  make_option(c("--data_type"),
              type = "character",
              default = NULL,
              help = "Data type used to generate qpAdm analysis",
              metavar = "character"),
  make_option(c("--model_name"),
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
arguments = parse_args(opt_parser, positional_arguments=TRUE);
opt <- arguments$options

if (is.null(opt$in_dir)){
  print_help(opt_parser)
  stop("In DIR for branch average FST RDS files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}


#' Creating a list of the rds files for all the replicates for each target
fst_reps_SummaryTable = data.table()
in_files_list = list.files(opt$in_dir)
print(in_files_list)
for(i in 1:length(in_files_list)){
  fst_rep = data.table(fread(file = paste0(opt$in_dir, in_files_list[i])))
  #fst_rep = data.table(fread(file = paste0(in_dir, in_files_list[i]))) ##! TESTING
  fst_reps_SummaryTable = rbind(fst_reps_SummaryTable, fst_rep)
}

fst_reps_SummaryTable$data_type <- opt$data_type
fst_reps_SummaryTable$model_name <- opt$model_name

# Save data table as RDS
#saveRDS(fst_reps_SummaryTable, file=opt$out)
write.table(fst_reps_SummaryTable, file = opt$out, sep = "\t", row.names = FALSE)


#' ####################################################################################################################################################################################################
#' TESTING
# setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/project__complex_msp_simulation/")
# in_dir = "out/ROAR_Output/simple_Demography_5kSims/branch_FST/"
# in_files_list = list.files(in_dir)
# fst_rep = data.table(fread(file = paste0(in_dir, in_files_list[1])))
