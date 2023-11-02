#' Generate qpAdm Summary Table
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(admixtools)
library(tidyverse)
library(optparse)

# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("-d", "--in_dir"),
              type = "character",
              default = NULL,
              help = "In Dir for qpAdm RDS files",
              metavar = "character"),
  make_option(c("-t","--data_type"),
              type = "character",
              default = NULL,
              help = "Data type used to generate qpAdm analysis",
              metavar = "character"),
  make_option(c("-p","--target_pop"),
              type = "character",
              default = NULL,
              help = "qpAdm Target Pop",
              metavar = "character"),
  make_option(c("p","--processing_scr"),
              type = "character",
              default = NULL,
              help = "Rscript for processing the qpAdm tables",
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
arguments = parse_args(opt_parser, positional_arguments=TRUE);
opt <- arguments$options

if (is.null(opt$in_dir)){
  print_help(opt_parser)
  stop("In DIR for qpAdm RDS files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}

#' Parameters
modelName = opt$model_name # simulation model name
dataType= "simulated_F2_branch"
DIR_Path = opt$in_dir # input data DIR
target = opt$target_pop

message("TARGET POPULATION ::::: ", target)

#' Generate list of sim qpAdm RDS files (keeping alpha numeric order)
file_list <- list.files(path = DIR_Path, pattern = paste0(target, "\\.rds$"), full.names = TRUE)
sorted_index <- order(as.numeric(gsub("[^0-9]", "", file_list)))
file_list <- file_list[sorted_index]
inDIR = c()
for(f in file_list){
  inDIR[[length(inDIR) +1]] <- f
}

# Load admixtools2_qpAdm_rotation_processing.R script.
source(paste0(opt$processing_scr))

# Generate Summary Table
qpAdm_rotation_summary_table = qpAdm_rotate_processing(DIR = inDIR, sim_model = modelName, dataType=dataType)

# Write output Table
write.table(qpAdm_rotation_summary_table, file=opt$out, row.names = FALSE)
