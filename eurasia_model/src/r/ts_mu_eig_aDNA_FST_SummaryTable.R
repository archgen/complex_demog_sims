#' Generate FST sim aDNA Summary Table
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
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
args <- arguments$args

if (is.null(opt$in_dir)){
  print_help(opt_parser)
  stop("In DIR for admixture F3 RDS files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}

#' Creating a list of the rds files for all the replicates for each target
fst_reps_SummaryTable = data.table()
in_rds_list = list.files(opt$in_dir)
print(in_rds_list)
for(i in 1:length(in_rds_list)){
  fst_rep = data.table(readRDS(file = paste0(opt$in_dir, in_rds_list[i])))
  #fst_rep = data.table(readRDS(file = paste0(in_dir, in_rds_list[i]))) ##! TESTING
  fst_DT = fst_rep[!duplicated(fst_rep[, c("est", "se", "rep")]), ]
  fst_DT = data.table(fst_DT)[pop1 != pop2]
  fst_reps_SummaryTable = rbind(fst_reps_SummaryTable, fst_DT)
}

# Save data table as RDS
#saveRDS(fst_reps_SummaryTable, file=opt$out)
write.table(fst_reps_SummaryTable, file = opt$out, sep = "\t", row.names = FALSE)
#' ####################################################################################################################################################################################################


