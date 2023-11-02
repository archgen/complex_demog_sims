#' Compute admixtools2 analyses (qpAdm)
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
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
targetPop = p_list[grepl("T", p_list)]
sourcePops = p_list[grepl("S|R", p_list)]
out = opt$out

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
# Define qpAdm rotation models
qpAdm_Rotate_Models <- rotate_models(leftright=c(unique(sourcePops)), target=targetPop, rightfix = NULL)
qpAdm_Rotate_Models = qpAdm_Rotate_Models %>% filter(lengths(left) <= 2)
# Run qpAdm rotation
qpadm_multi_run=qpadm_multi(f2_branch_matrix, qpAdm_Rotate_Models, full_results=TRUE, fudge_twice=TRUE)
# Save output
saveRDS(qpadm_multi_run, file=out)
