#' Compute admixtools2 analyses (admixture F3) on simulated aDNA data
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
  make_option(c("--sim_geno"),
              type = "character",
              default = NULL,
              help = "geno simulated eigenstrat file",
              metavar = "character"),
  make_option(c("--samples"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
              metavar = "character"),
  make_option(c("--target"),
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
  make_option(c("--replicate"),
              type = "character",
              default = NULL,
              help = "Iteration or replicate number",
              metavar = "integer"),
  make_option(c("--out"),
              type = "character",
              default = NULL,
              help = "output FILENAME & PATH",
              metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$sim_geno)){
  print_help(opt_parser)
  stop("simulated genotype files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}


# Parameters
sampleSheet = data.table(read.csv(opt$samples))
sampleSheet[sampleSheet == ''] <- NA
targetPop = opt$target
targetPop = gsub("\\[|\\]", "", targetPop)
SourceRightPops = sampleSheet[sampleSheet[, qpAdmRotate %in% c("source", "right")]]$Pop
out = opt$out
eig_prefix = str_split(opt$sim_geno, ".geno")[[1]][1]
iter = opt$replicate

#' Run admixture F3
aDNA_admixf3 = qp3pop(eig_prefix, targetPop, SourceRightPops, SourceRightPops, blgsize = 4e6)
aDNA_admixf3$rep <- iter # add simulation replicate information
aDNA_admixf3$data_type <- opt$data_type # add simulation data type information
aDNA_admixf3$model_name <- opt$model_name # add simulation name information


#' Save output
saveRDS(aDNA_admixf3, file=out)

