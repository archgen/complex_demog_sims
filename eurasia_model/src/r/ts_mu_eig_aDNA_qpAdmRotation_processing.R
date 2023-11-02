#' Process mutated tree-sequence with ancient DNA conditions
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(vroom)
library(tidyverse)
library(optparse)
library(tictoc)
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
  make_option(c("--real_eig_prefix"),
              type = "character",
              default = NULL,
              help = "Real aDNA eigenstrat prefix & directory",
              metavar = "character"),
  make_option(c("--sample_sheet"),
              type = "character",
              default = NULL,
              help = "Input sample sheet",
              metavar = "character"),
  make_option(c("--recom_map"),
              type = "character",
              default = NULL,
              help = "Input recombination_map",
              metavar = "character"),
  make_option(c("--script_func"),
              type = "character",
              default = NULL,
              help = "aDNA conditions functions R script",
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
  stop("In DIR with mutated tree-sequences must be provided (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out FILENAME & PATH must be supplied (output file).n", call. = FALSE)
}

# ### TESTING
# sampleInfo <- fread("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/input/sampleSheet_HomSap__AncientEurasia_9K19.txt", fill=TRUE)
# RscriptDIR <- c("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/src/R/")
# outData_pseudohap_fn <- c("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/out/AncientEurasia_9K19/ts_mu_eig/ts_mu_AncientEurasia_9K19_rep_1_ascertain_pseudohap_1240K_missing")
# simData_fnprefix <- c("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/out/AncientEurasia_9K19/ts_mu_eig/ts_mu_AncientEurasia_9K19_rep_1")
# realData_fnprefix <- c("1240K_public_v50.0_SimsEurasia")
# realData_dir <- c("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/input/")
# source("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/src/R/Processing_aDNAConditions.R")
# recom_map <- c("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/input/hapmap/recom_map_HapMapII_GRCh37__Chrs21_22.txt")
#setwd("/Users/mkw5910/Documents/PSU_sims/snakemake/project__complex_msp_simulation/")
#sim_eigenstrat = read_eigenstrat(eig_geno="out/HomSap__AncientEurasia_9K19/ts_mu_eig/ts_mu_HomSap__AncientEurasia_9K19_rep_5.geno", 
#                                 eig_snp="out/HomSap__AncientEurasia_9K19/ts_mu_eig/ts_mu_HomSap__AncientEurasia_9K19_rep_5.snp", 
#                                 eig_ind="out/HomSap__AncientEurasia_9K19/ts_mu_eig/ts_mu_HomSap__AncientEurasia_9K19_rep_5.ind") 
# ### TESTING



################################################
# Eigenstrat aDNA conditions Processing        #
################################################
# Load eigenstrat_aDNA_conditions_processing.R script.
source(paste0(opt$script_func))

# Load simulated & real data
message("Loading simulated data")
tic(":::: Loading simulated data Run Time = ")
sim_eigenstrat = read_eigenstrat(eig_geno=(opt$sim_geno), eig_snp=(opt$sim_snp), eig_ind=(opt$sim_ind))  # Large file, this takes ~ 1-2min
toc()
message("Loading real data")
tic(":::: Loading real data Run Time = ")
real_eigenstrat = read_eigenstrat(eig_geno=paste0(opt$real_eig_prefix, ".geno"),
                                                  eig_snp=paste0(opt$real_eig_prefix, ".snp"),
                                                  eig_ind=paste0(opt$real_eig_prefix, ".ind")) # Data must be in Eigenstrat format not AncestryMap

toc()
message("Loading hapmap data")
tic(":::: Loading hapmap Run Time = ")
recom_map <- vroom(opt$recom_map, show_col_types = FALSE)
toc()


# Define ancient and modern populations
sampleInfo <- fread(opt$sample_sheet, fill=TRUE)
eig_sampleInfo = copy(sampleInfo %>% filter(Output == "eig"))
ancientPops <- c(eig_sampleInfo[eig_sampleInfo$Generations != 0, Pop])
modernPops <- c(eig_sampleInfo[eig_sampleInfo$Generations == 0, Pop])
ascertainInds <- eig_sampleInfo[grepl("ascertain", eig_sampleInfo$Pop, fixed = TRUE), Pop]


# Define the individuals to ascertain on.
ascertainment_Inds = lapply(c(ascertainInds), function(x) which(sim_eigenstrat$ind$group == x))
#' Print to terminal the list of ascertainment individuals and their position in the genotype files
cat("Ascertainment Inds & positions in .geno file: ")
for (i in seq_along(ascertainInds)) {
  cat(paste(ascertainInds[[i]], " ::: ", ascertainment_Inds[[i]], "   ", sep = ""))
}
cat("\n")
#' Print to terminal the list of ancient and present-day populations
cat(paste((" ** Ancient populations ** ")), paste0("\n", ancientPops, collapse = ""))
cat(paste((" ** Present-day populations ** ")), paste0("\n", modernPops, collapse = ""))


# Run functions on data
"Running Functions ... "
message("Filter tri-allelic sites")
tic(":::: Tri-allelic filtering Run Time = ")
sim_eigenstrat_biAllele = triAllelFilterd(eigenstrat_file = sim_eigenstrat)
toc()
message("Adding ascertaiment")
tic(":::: Ascertainment Run Time = ")
sim_eigenstrat_ascertained = add_ascertainment_bias(eigenstrat_file = sim_eigenstrat_biAllele, ascertainment_based_on = ascertainment_Inds)
toc()
message("Making diploid")
tic(":::: Diploidy Run Time = ")
sim_eigenstrat_ascertained_diploid = makeSimDiploid(sim_eigenstrat_ascertained)
toc()
who_gets_pseudoHap = which(sim_eigenstrat_ascertained_diploid$ind$group %in% ancientPops)
message("Pseudohaploidising data")
tic(":::: PseudoHap Conversion Run Time = ")
sim_eigenstrat_ascertained_pseudohaploid = pseudoHaplodize(eigenstrat_file = sim_eigenstrat_ascertained_diploid, who_gets_pseudoHap = who_gets_pseudoHap)
toc()
who_gets_missing = which(sim_eigenstrat_ascertained_pseudohaploid$ind$group %in% ancientPops)
message("Downsampling to 1240K & adding missing data")
tic(":::: Downsampling & Missing Data Process Run Time = ")
sim_eigenstrat_ascertained_pseudohaploid_1240K_missing = v50_0_1240Kdownsamp_add_missing(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid, reference_file = real_eigenstrat, who_gets_missing = who_gets_missing)
toc()
message("recode .snp file")
tic(":::: Recoding snp file Run Time = ")
sim_eigenstrat_ascertained_pseudohaploid_1240K_missing_recodeSNP = recode_snpFile(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid_1240K_missing, hapmap_file = recom_map)
#!! >> When you get the new code on cum_map_MSP add this to the source script
toc()
message("Writing pseudohaploid eig data ...")
tic(":::: Writing pseudoHap Data Run Time = ")
write_eigenstrat_pshap(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid_1240K_missing_recodeSNP, fn = opt$out)
toc()



