#' Compute admixtools2 analyses (qpAdm Rotation on ts F2 branch lengths)
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#' Set path to R library dir
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
  make_option(c("-s", "--samples"),
              type = "character",
              default = NULL,
              help = "Sample sheet",
              metavar = "character"),
  make_option(c("-t", "--target"),
              type = "character",
              default = NULL,
              help = "qpAdm Target population",
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

# TESTING
#sampleSheet = data.table(read.csv("../Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/input/Snakemake_sampleSheet_ts_aDNA_complexDemography.txt"))

# Parameters
f2_branch_matrix = readRDS(opt$f2_matrix, refhook = NULL)
sampleSheet = data.table(read.csv(opt$samples))
sampleSheet[sampleSheet == ''] <- NA
targetPop = opt$target
targetPop = gsub("\\[|\\]", "", targetPop)
sourcePops = sampleSheet %>% filter(qpAdmRotate == "source")
rightPops = sampleSheet %>% filter(qpAdmRotate == "right")
out = opt$out


# qpAdm Analysis

# Define qpAdm rotation models
right_pops = unique(rightPops$Pop)
right_pops = right_pops[!grepl("Mbuti", right_pops)]
qpAdm_Rotate_Models <- rotate_models(leftright=c(unique(sourcePops$Pop)),
                                     target=targetPop,
                                     rightfix = right_pops)
qpAdm_Rotate_Models = qpAdm_Rotate_Models %>% filter(lengths(left) <= 3)

left_sets = do.call(rbind, lapply(qpAdm_Rotate_Models[[1]], function(x) {c(x, rep(NA, max(sapply(qpAdm_Rotate_Models[[1]], length)) - length(x)))}))
right_sets <- do.call(rbind, lapply(qpAdm_Rotate_Models[[2]], function(x) {c(x, rep(NA, max(sapply(qpAdm_Rotate_Models[[2]], length)) - length(x)))}))

# Check if the length of left_sets and right_sets is equal
if (nrow(left_sets) != nrow(right_sets)) {
  warning("Length of left and right sets is unequal")
  # Add any additional actions you want to take if the lengths are unequal
} else {
  # Run qpAdm rotation
  qpadm_rotate_run = list()
  qpadm_out_1 = for(i in 1:nrow(left_sets)) {
    pops_left = c(left_sets[i,][!is.na(left_sets[i,])])
    pops_right = c("Mbuti", right_sets[i,][!is.na(right_sets[i,])])
    #message("Left pops = ", "model ", i, "   ", paste0(left_pops, collapse = ", "))
    #message("Right pops = ", "model ", i, "   ", paste0(right_pops, collapse = ", "))

    out_qpadm = qpadm(data = f2_branch_matrix,
                target = targetPop, left = pops_left,
                right = pops_right,
                boot=F, fudge_twice=T, getcov=T, cpp=T, verbose=T)
    #nm = paste0(targetPop, "; ", paste(left_pops, collapse=", "), "; ", paste(right_pops, collapse=", "))
    qpadm_rotate_run[[i]] = out_qpadm
  }
}


# Save output
saveRDS(qpadm_rotate_run, file=out)


## qpAdm Analysis
#message("")
#message("Beginning qpAdm analysis on target population ", targetPop)
## Define qpAdm rotation models
#qpAdm_Rotate_Models <- rotateMax3Left(leftright=c(unique(sourcePops$Pop)), target=targetPop, rightfix = NULL)
#qpAdm_Rotate_Models <- rotate_models(leftright=c(unique(sourcePops$Pop)), target=targetPop, rightfix = c(unique(rightPops$Pop)))
#qpAdm_Rotate_Models_filt <- qpAdm_Rotate_Models %>% filter(lengths(left) <= 2)
## Run qpAdm rotation
#qpadm_multi_run=qpadm_multi(f2_branch_matrix, qpAdm_Rotate_Models_filt, full_results=TRUE, fudge_twice = TRUE)
## Save output
#saveRDS(qpadm_multi_run, file=out)
