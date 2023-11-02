#' Generate a matrix of the tree_sequence f2 blocks
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#.libPaths("/storage/home/m/mkw5910/.conda/envs/snakemake_base/lib/R/library")
library(tidyverse)
library(data.table)
library(optparse)

# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("-t", "--f2_table"),
              type = "character",
              default = NULL,
              help = "msprime tree sequence path",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "output file name and path",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$f2_table)){
  print_help(opt_parser)
  stop("Tree sequence must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out directory and filename must be supplied (output file).n", call. = FALSE)
}

f2_table = fread(opt$f2_table)
names <- unique(unlist(strsplit(colnames(f2_table), "_AND_")))
num_windows = nrow(f2_table)
ncol <- length(names)
nrow <- length(names)
dim=c(ncol, nrow, num_windows)
column.names <- names
row.names <- names
matrix.names <- rep(paste0("block", 1:num_windows), 1)
simF2_blocks <- array(0, dim=dim, dimnames=list(row.names, column.names, matrix.names))

#' Populate matrix with values from data table 
N <- strsplit(colnames(f2_table), "_AND_")
for(i in 1:length(N)){
  for (k in 1:num_windows) {
    p1 = N[i][[1]][[1]]
    p2 = N[i][[1]][[2]]
    simF2_blocks[p2, p1, k] <- as.numeric(f2_table[k, ..i])
    simF2_blocks[p1, p2, k] <- as.numeric(f2_table[k, ..i])
  }
}

#' Write f2 matrix out
write_rds(simF2_blocks, opt$out)
