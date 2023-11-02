#' Matthew Williams: 
#'  
#' Load in libraries
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(admixtools)
library(tidyverse)
library(optparse)
library(data.table)

# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("-m", "--f2_matrix"),
              type = "character",
              default = NULL,
              help = "F2 branch matrix",
              metavar = "character"),
  make_option(c("-r", "--replicate"),
              type = "character",
              default = NULL,
              help = "Iteration or replicate number",
              metavar = "integer"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "output FILENAME & PATH",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$f2_matrix)){
  print_help(opt_parser)
  stop("F2 branch matrix must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}


# Parameters
f2_branch_matrix = readRDS(opt$f2_matrix, refhook = NULL)
p_list = colnames(f2_branch_matrix)
targetPop =  p_list[grepl("T", p_list)]
SourceRightPops = p_list[grepl("S|R", p_list)]
out = opt$out
iter = opt$replicate

#' Run admixture F3
branch_admixf3 = qp3pop(f2_branch_matrix, targetPop, SourceRightPops, SourceRightPops)
branch_admixf3$rep <- iter # add simulation replicate information
#' Save output
saveRDS(branch_admixf3, file=out)