#' Compute admixtools2 qpAdm Rotation on simulated aDNA data
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
  make_option(c("--sim_snp"),
              type = "character",
              default = NULL,
              help = "snp simulated eigenstrat file",
              metavar = "character"),
  make_option(c("--sim_ind"),
              type = "character",
              default = NULL,
              help = "ind simulated eigenstrat file",
              metavar = "character"),
  make_option(c("--sample_sheet"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
              metavar = "character"),
  make_option(c("-t", "--target"),
              type = "character",
              default = NULL,
              help = "Target pop",
              metavar = "character"),
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
targetPop = opt$target
targetPop = gsub("\\[|\\]", "", targetPop)
sourcePops = sampleSheet %>% filter(qpAdmRotate == "source")
rightPops = sampleSheet %>% filter(qpAdmRotate == "right")
out = opt$out
eig_prefix = str_split(opt$sim_geno, ".geno")[[1]][1]


##################################################################
# Functions to create qpAdm Models with Max 3 source populations #
##################################################################
all_lr2 = function(pops, rightfix = 0) {
  # splits pops into all possible combinations of left and right
  maxleft = min(3, (length(pops) -1) /2 )
  left = map(1:maxleft, ~combn(pops, ., simplify = F)) %>% flatten
  right = map(left, ~setdiff(pops, .))
  namedList(left, right)
}
# min(length(pops), floor((length(pops) + rightfix - 1)/2))
rotateMax3Left = function (leftright, target, rightfix = NULL){
  all_lr2(leftright, length(rightfix)) %>% as_tibble %>% rowwise %>%
    mutate(right = list(c(right, rightfix)), target = target) %>%
    ungroup
}

# qpAdm Analysis
message("")
message("Beginning qpAdm analysis on target population ", targetPop)
message("Using Source populations ", unique(sourcePops$Pop))
message("And Right populations ", unique(rightPops$Pop))
message("on the dataset  ", eig_prefix)
message("Writing the output to   ", out)


# Define qpAdm rotation models
#qpAdm_Rotate_Models <- rotateMax3Left(leftright=c(unique(sourcePops$Pop)), target=targetPop, rightfix = c(unique(rightPops$Pop)))
qpAdm_Rotate_Models <- rotate_models(leftright=c(unique(sourcePops$Pop)), target=targetPop, rightfix = c(unique(rightPops$Pop)))
qpAdm_Rotate_Models_filt <- qpAdm_Rotate_Models %>% filter(lengths(left) <= 3)
# Run qpAdm rotation
qpadm_multi_run=qpadm_multi(eig_prefix, qpAdm_Rotate_Models_filt, allsnps = TRUE, full_results=TRUE, fudge_twice = TRUE)
# Save output
saveRDS(qpadm_multi_run, file=out)
