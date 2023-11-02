#' Compute admixtools2 analyses (FST) on simulated aDNA data
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
  make_option(c("--sample_sheet"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
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
sampleSheet = data.table(read.csv(opt$sample_sheet))
sampleSheet[sampleSheet == ''] <- NA
TargetSourceRightPops = sampleSheet[sampleSheet[, qpAdmRotate %in% c("source", "right", "target")]]$Pop
out = opt$out
eig_prefix = str_split(opt$sim_geno, ".geno")[[1]][1]
iter = opt$replicate

#' Run FST
aDNA_fst = fst(eig_prefix, c(TargetSourceRightPops), c(TargetSourceRightPops),
               all_snps = TRUE, adjust_pseudohaploid = TRUE, fudge_twice=TRUE, blgsize = 4e6)
aDNA_fst$rep <- iter # add simulation replicate information
aDNA_fst$data_type <- opt$data_type # add simulation data type information
aDNA_fst$model_name <- opt$model_name # add simulation model name information


#' Save output
saveRDS(aDNA_fst, file=out)
