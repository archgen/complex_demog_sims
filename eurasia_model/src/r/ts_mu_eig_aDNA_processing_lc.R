#' Process mutated tree-sequence with ancient DNA conditions
#' 
#' Matthew Williams: 
#'  

####################################
# PACKAGES                         #
####################################
# Load the dev data.table package from the specified library location
library(data.table, lib.loc="./dev_RLibs")
# Load the vroom package from the specified library location
#library(vroom, lib.loc = "/storage/home/mkw5910/.conda/envs/snakemake_base/lib/R/library")
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(vroom)
library(tidyverse)
library(optparse)
library(tictoc)

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
  make_option(c("--aDNA_data_dir"),
              type = "character",
              default = NULL,
              help = "Empirical aDNA directory",
              metavar = "character"),
  make_option(c("--empir_data_prefix"),
              type = "character",
              default = NULL,
              help = "Empirical aDNA data prefix",
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
  make_option(c("--length"),
              type = "character",
              default = NULL,
              help = "Simulated genome length (chr22 or WG)",
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


################################################
# Eigenstrat aDNA conditions Processing        #
################################################
# Load eigenstrat_aDNA_conditions_processing.R script.
source(paste0(opt$script_func))

# Load simulated & defining real data input
message("Loading simulated data", rep("\n", 1), sep = "")
tic(":::: Loading simulated data Run Time = ")
sim_eigenstrat = read_eigenstrat(eig_geno=(opt$sim_geno), eig_snp=(opt$sim_snp), eig_ind=(opt$sim_ind))  # Large file, this takes ~ 1-2min
# Input (data dir + empir data prefix)
empir_eig_input = paste0(opt$aDNA_data_dir, "/", opt$empir_data_prefix)
message("", rep("\n", 10), sep = "")
message("EMPERICAL DATA DIRECTORY")
print(empir_eig_input)
message("", rep("\n", 10), sep = "")
toc()
message("", rep("\n", 10), sep = "")

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
message("", rep("\n", 10), sep = "")
#' Print to terminal the list of ancient and present-day populations
cat(paste((" ** Ancient populations ** ")), paste0("\n", ancientPops, collapse = ""))
message("", rep("\n", 3), sep = "")
cat(paste((" ** Present-day populations ** ")), paste0("\n", modernPops, collapse = ""))
message("", rep("\n", 10), sep = "")


# Run functions on data
message("**** Running Processing Functions ****  ", rep("\n", 3), sep = "")

message("Filter tri-allelic sites", rep("\n", 0), sep = "")
tic(":::: Tri-allelic filtering Run Time = ")
sim_eigenstrat_biAllele = triAllelFilterd(eigenstrat_file = sim_eigenstrat)
rm(sim_eigenstrat)
toc()
message("", rep("\n", 3), sep = "")


message("Adding ascertaiment", rep("\n", 0), sep = "")
tic(":::: Ascertainment Run Time = ")
sim_eigenstrat_ascertained = add_ascertainment_bias(eigenstrat_file = sim_eigenstrat_biAllele, ascertainment_based_on = ascertainment_Inds)
rm(sim_eigenstrat_biAllele)
toc()
message("", rep("\n", 3), sep = "")


message("Making diploid", rep("\n", 0), sep = "")
tic(":::: Diploidy Run Time = ")
sim_eigenstrat_ascertained_diploid = makeSimDiploid(sim_eigenstrat_ascertained)
rm(sim_eigenstrat_ascertained)
toc()
message("", rep("\n", 3), sep = "")


message("Pseudohaploidising data", rep("\n", 0), sep = "")
tic(":::: PseudoHap Conversion Run Time = ")
who_gets_pseudoHap = which(sim_eigenstrat_ascertained_diploid$ind$group %in% ancientPops)
sim_eigenstrat_ascertained_pseudohaploid = pseudoHaplodize(eigenstrat_file = sim_eigenstrat_ascertained_diploid, who_gets_pseudoHap = who_gets_pseudoHap)
rm(sim_eigenstrat_ascertained_diploid)
toc()
message("", rep("\n", 3), sep = "")


message("Downsampling to 1240K & adding missing data", rep("\n", 1), sep = "")
tic(":::: Downsampling & Missing Data Process Run Time = ")
message("Loading real data")
tic(":::: Loading real data Run Time = ")
real_eigenstrat = read_eigenstrat(eig_geno=paste0(empir_eig_input, ".geno"),
                                                  eig_snp=paste0(empir_eig_input, ".snp"),
                                                  eig_ind=paste0(empir_eig_input, ".ind")) # Data must be in Eigenstrat format not AncestryMap

toc()
message("", rep("\n", 3), sep = "")

who_gets_missing = which(sim_eigenstrat_ascertained_pseudohaploid$ind$group %in% ancientPops)
sim_eigenstrat_ascertained_pseudohaploid_1240K_missing = v50_0_1240Kdownsamp_add_missing(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid, reference_file = real_eigenstrat, who_gets_missing = who_gets_missing)
rm(sim_eigenstrat_ascertained_pseudohaploid)
rm(real_eigenstrat)
toc()
message("", rep("\n", 3), sep = "")

length=opt$length
if(length == "WG") {
    message("recode .snp file")
    tic(":::: Recoding snp file Run Time = ")
    sim_eigenstrat_ascertained_pseudohaploid_1240K_missing_recodeSNP = recode_snpFile(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid_1240K_missing)
    toc()
    rm(sim_eigenstrat_ascertained_pseudohaploid_1240K_missing)
} else if(length == "chr22") {
    sim_eigenstrat_ascertained_pseudohaploid_1240K_missing_recodeSNP = sim_eigenstrat_ascertained_pseudohaploid_1240K_missing
}


message("Writing pseudohaploid eig data ..", rep("\n", 2), sep = "")
tic(":::: Writing pseudoHap Data Run Time = ")
write_eigenstrat_pshap(eigenstrat_file = sim_eigenstrat_ascertained_pseudohaploid_1240K_missing_recodeSNP, fn = opt$out)
toc()
message("", rep("\n", 10), sep = "")