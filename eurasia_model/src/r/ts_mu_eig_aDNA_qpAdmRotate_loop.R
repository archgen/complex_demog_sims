#' Compute admixtools2 qpAdm Rotation on simulated aDNA data
#' 
#' Matthew Williams: 
#'  
#' Load in libraries

#### TESTING
#sim_geno = "/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/out/AncientEurasia_9K19/ts_mu_eig_aDNA/ts_mu_AncientEurasia_9K19_rep_1_ascertain_pseudohap_1240K_missing.geno"
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
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
outfile = opt$out
eig_prefix = str_split(opt$sim_geno, ".geno")[[1]][1]

#print(eig_prefix)
#ind = fread(opt$sim_ind, header = FALSE)
#ind_prefix = paste0(eig_prefix, ".ind")
#ind=fread(ind_prefix, header = FALSE) 
#print(ind)

# qpAdm Analysis
message("")
message("Beginning qpAdm analysis on target population ", targetPop)
message("Using Source populations ", unique(sourcePops$Pop))
message("And Right populations ", unique(rightPops$Pop))
message("on the dataset  ", eig_prefix)
message("Writing the output to   ", outfile)


# Define qpAdm rotation models
right_pops = unique(rightPops$Pop)
right_pops = right_pops[!grepl("Mbuti", right_pops)]
qpAdm_Rotate_Models = rotate_models(leftright=c(unique(sourcePops$Pop)), 
                                     target=targetPop, 
                                     rightfix = right_pops)
qpAdm_Rotate_Models = qpAdm_Rotate_Models %>% filter(lengths(left) <= 3)

left_sets = do.call(rbind, lapply(qpAdm_Rotate_Models[[1]], function(x) {c(x, rep(NA, max(sapply(qpAdm_Rotate_Models[[1]], length)) - length(x)))}))
right_sets = do.call(rbind, lapply(qpAdm_Rotate_Models[[2]], function(x) {c(x, rep(NA, max(sapply(qpAdm_Rotate_Models[[2]], length)) - length(x)))}))

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
    
    out = qpadm(data = eig_prefix, 
                target = targetPop, left = pops_left, 
                right = pops_right, 
                boot=F, fudge_twice=T, getcov=T, cpp=T, verbose=T, allsnps=T, blgsize = 4e6)
    #nm = paste0(target_pop, "; ", paste(left_pops, collapse=", "), "; ", paste(right_pops, collapse=", "))
    qpadm_rotate_run[[i]] = out
  }
}


# Save output
saveRDS(qpadm_rotate_run, file=outfile)
