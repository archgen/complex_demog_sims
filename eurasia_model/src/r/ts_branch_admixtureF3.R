#' Compute admixtools2 analyses (admixture f3)
#' 
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
  make_option(c("-s", "--samples"),
              type = "character",
              default = NULL,
              help = "Sample sheet",
              metavar = "character"),
  make_option(c("-t", "--target"),
              type = "character",
              default = NULL,
              help = "Target pop",
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
sampleSheet = data.table(read.csv(opt$samples))
sampleSheet[sampleSheet == ''] <- NA
targetPop = opt$target
targetPop = gsub("\\[|\\]", "", targetPop)
SourceRightPops = sampleSheet[sampleSheet[, qpAdmRotate %in% c("source", "right")]]$Pop
out = opt$out
iter = opt$replicate

#' Run admixture F3
branch_admixf3 = qp3pop(f2_branch_matrix, targetPop, SourceRightPops, SourceRightPops)
branch_admixf3$rep <- iter # add simulation replicate information
branch_admixf3$data_type <- opt$data_type # add simulation data type information
branch_admixf3$model_name <- opt$model_name # add simulation name information

#' Save output
saveRDS(branch_admixf3, file=out)
