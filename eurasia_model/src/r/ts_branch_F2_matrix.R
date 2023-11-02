#' Generate a matrix of the tree_sequence f2 blocks
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
#.libPaths("/storage/group/cdh5313/default/mkw5910/sims/complex_Demography/.snakemake/conda/5af5b266/lib/R/library/") # conda env R lib
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

#' Testing
# f2_table = fread("./out/ROAR_Output/simple_Demography_5kSims/TEST_F2_array.txt")
# f2_table = data.table("p1_AND_p2" = c(0.1, 0.2, 0.3), 
#                  "p1_AND_p3" = c(0.4, 0.5, 0.6),
#                  "p2_AND_p3" = c(0.9, 0.8, 0.7))



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
