#' Compute admixtools2 analyses (qpAdm Replacement simulated aDNA data)
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#' Set path to R library dir
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")

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
              help = "F2 branch matrix",
              metavar = "character"),
  make_option(c("--samples"),
              type = "character",
              default = NULL,
              help = "Sample sheet",
              metavar = "character"),
  make_option(c("--iteration"),
              type = "numeric",
              default = NULL,
              help = "Sim Replicate",
              metavar = "numeric"),
  make_option(c("-t", "--target"),
              type = "character",
              default = NULL,
              help = "qpAdm Target population",
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
  stop("Simulated genotype files must be supplied (input file).n", call. = FALSE)
}
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Out list must be supplied (output file).n", call. = FALSE)
}


#' TESTING TESTING 
#setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/out/complexDemography_model_A_HapMapII_GRCh37__Chrs21_22/")
#eig_prefix = "./ts_mu_eig_aDNA/ts_mu_complexDemography_model_A_HapMapII_GRCh37__Chrs21_22_rep_1_ascertain_pseudohap_1240K_missing"
#f2_branch_matrix <- readRDS("./complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource/admixtools2/ts_branch_F2_matrix/ts_branch_F2_matrix__complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource_rep_1.rds")
#sampleSheet <- data.table(read_csv("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/project__complex_msp_simulation/input/sampleSheet_complexSimsModel__AncientEurasia_A.txt"))
#' TESTING multiple target pops target_pops = c(target_pops, "ascertain_CHB__gen_0", "ascertain_AFR__gen_0")
#ind = fread("./complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource/ts_mu_eig_aDNA/ts_mu_complexSimsModel__AncientEurasia_A_HapMapII_GRCh37__Chrs21_22_wMediSource_rep_1_ascertain_pseudohap_1240K_missing.ind",
#            header=FALSE)

# Parameters
rep = opt$iteration
eig_prefix = str_split(opt$sim_geno, ".geno")[[1]][1]
sampleSheet = data.table(read.csv(opt$samples))
sampleSheet[sampleSheet == ''] <- NA
targetPop = opt$target
targetPop = gsub("\\[|\\]", "", targetPop)
sourcePops = sampleSheet %>% filter(qpAdmFixed == "source")
rightPops = sampleSheet %>% filter(qpAdmFixed == "right")


#' Running qpAdm Replacement for three & two and single source models 
qpAdm_Fixed_sources_summary <- data.frame("target" = as.character(),
                                          "source.p1"= as.character(), "source.p2"= as.character(),  "source.p3"= as.character(), 
                                          "weight.p1"= as.numeric(), "weight.p2"= as.numeric(),"weight.p3"= as.numeric(),
                                          "p.value"= as.numeric(), 
                                          "se.p1"= as.numeric(), "se.p2"= as.numeric(), "se.p3"= as.numeric(),  
                                          "f4rank"= as.numeric(),
                                          "rep" = as.numeric()) 

if(length(rightPops$Pop) > 3){
  message("Running 3-2-1 Source qpAdm models")
  #' Three-way models
  #' Generate all unique combinations of 3 elements
  pop3waycombinations <- combn(sourcePops$Pop, 3)
  # Convert the matrix into a list of lists
  pop3waycom_list <- lapply(seq_len(ncol(pop3waycombinations)), function(i) pop3waycombinations[, i])
  
  for(i in 1:length(pop3waycom_list)){
    message("For target populations : ", targetPop, 
            "  Running 3 Source pairs  ", pop3waycom_list[[i]][1], " ", pop3waycom_list[[i]][2], " ", pop3waycom_list[[i]][3]) 
    qpAdm_result <- qpadm(eig_prefix, c(pop3waycom_list[[i]][1], pop3waycom_list[[i]][2], pop3waycom_list[[i]][3]), 
                          unique(c("Mbuti", c(rightPops$Pop))), targetPop, fudge_twice = TRUE, blgsize = 4e6) 
    qpAdm_Fixed_sources_summary[nrow(qpAdm_Fixed_sources_summary)+1, ] <- c(as.character(qpAdm_result$weights$target[1]),
                                                                            as.character(qpAdm_result$weights$left[1]),
                                                                            as.character(qpAdm_result$weights$left[2]),
                                                                            as.character(qpAdm_result$weights$left[3]),
                                                                            as.numeric(qpAdm_result$weights$weight[1]),
                                                                            as.numeric(qpAdm_result$weights$weight[2]),
                                                                            as.numeric(qpAdm_result$weights$weight[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                            as.numeric(qpAdm_result$weights$se[1]),
                                                                            as.numeric(qpAdm_result$weights$se[2]),
                                                                            as.numeric(qpAdm_result$weights$se[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$f4rank[1]),
                                                                            as.numeric(rep))
  }
  
  #' Two-way models
  #' Generate all unique combinations of 3 elements
  pop2waycombinations <- combn(sourcePops$Pop, 2)
  # Convert the matrix into a list of lists
  pop2waycom_list <- lapply(seq_len(ncol(pop2waycombinations)), function(i) pop2waycombinations[, i])
  
  for(i in 1:length(pop2waycom_list)){
    message("For target populations : ", targetPop, 
            "  Running 2 Source pairs  ", pop2waycom_list[[i]][1], " ", pop2waycom_list[[i]][2]) 
    qpAdm_result <- qpadm(eig_prefix, c(pop2waycom_list[[i]][1], pop2waycom_list[[i]][2]), 
                          unique(c("Mbuti", c(rightPops$Pop))), targetPop, fudge_twice = TRUE, blgsize = 4e6) 
    qpAdm_Fixed_sources_summary[nrow(qpAdm_Fixed_sources_summary)+1, ] <- c(as.character(qpAdm_result$weights$target[1]),
                                                                            as.character(qpAdm_result$weights$left[1]),
                                                                            as.character(qpAdm_result$weights$left[2]),
                                                                            as.character(qpAdm_result$weights$left[3]),
                                                                            as.numeric(qpAdm_result$weights$weight[1]),
                                                                            as.numeric(qpAdm_result$weights$weight[2]),
                                                                            as.numeric(qpAdm_result$weights$weight[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                            as.numeric(qpAdm_result$weights$se[1]),
                                                                            as.numeric(qpAdm_result$weights$se[2]),
                                                                            as.numeric(qpAdm_result$weights$se[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$f4rank[1]),
                                                                            as.numeric(rep))
  }
  
  #' Single-source models
  for(i in (sourcePops$Pop)){
    message("For target populations : ", targetPop, 
            "  Running Single Source ", i) 
    qpAdm_result <- qpadm(eig_prefix, c(i), 
                          unique(c("Mbuti", c(rightPops$Pop))), targetPop, fudge_twice = TRUE, blgsize = 4e6) 
    qpAdm_Fixed_sources_summary[nrow(qpAdm_Fixed_sources_summary)+1, ] <- c(as.character(qpAdm_result$weights$target[1]),
                                                                            as.character(qpAdm_result$weights$left[1]),
                                                                            as.character(qpAdm_result$weights$left[2]),
                                                                            as.character(qpAdm_result$weights$left[3]),
                                                                            as.numeric(qpAdm_result$weights$weight[1]),
                                                                            as.numeric(qpAdm_result$weights$weight[2]),
                                                                            as.numeric(qpAdm_result$weights$weight[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                            as.numeric(qpAdm_result$weights$se[1]),
                                                                            as.numeric(qpAdm_result$weights$se[2]),
                                                                            as.numeric(qpAdm_result$weights$se[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$f4rank[1]),
                                                                            as.numeric(rep))
  }
  
  #' Return qpAdm summaries  
  #return(qpAdm_Fixed_sources_summary)
  
  
} else {
  message("Running 2-1 Source qpAdm models")
  #' Two-way models
  #' Generate all unique combinations of 3 elements
  pop2waycombinations <- combn(sourcePops$Pop, 2)
  # Convert the matrix into a list of lists
  pop2waycom_list <- lapply(seq_len(ncol(pop2waycombinations)), function(i) pop2waycombinations[, i])
  
  for(i in 1:length(pop2waycom_list)){
    message("For target populations : ", targetPop, 
            "  Running 2 Source pairs  ", pop2waycom_list[[i]][1], " ", pop2waycom_list[[i]][2]) 
    qpAdm_result <- qpadm(eig_prefix, c(pop2waycom_list[[i]][1], pop2waycom_list[[i]][2]), 
                          unique(c("Mbuti", c(rightPops$Pop))), targetPop, fudge_twice = TRUE, blgsize = 4e6) 
    qpAdm_Fixed_sources_summary[nrow(qpAdm_Fixed_sources_summary)+1, ] <- c(as.character(qpAdm_result$weights$target[1]),
                                                                            as.character(qpAdm_result$weights$left[1]),
                                                                            as.character(qpAdm_result$weights$left[2]),
                                                                            as.character(qpAdm_result$weights$left[3]),
                                                                            as.numeric(qpAdm_result$weights$weight[1]),
                                                                            as.numeric(qpAdm_result$weights$weight[2]),
                                                                            as.numeric(qpAdm_result$weights$weight[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                            as.numeric(qpAdm_result$weights$se[1]),
                                                                            as.numeric(qpAdm_result$weights$se[2]),
                                                                            as.numeric(qpAdm_result$weights$se[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$f4rank[1]),
                                                                            as.numeric(rep))
  }
  
  #' Single-source models
  for(i in (sourcePops$Pop)){
    message("For target populations : ", targetPop, 
            "  Running Single Source ", i) 
    qpAdm_result <- qpadm(eig_prefix, c(i), 
                          unique(c("Mbuti", c(rightPops$Pop))), targetPop, fudge_twice = TRUE, blgsize = 4e6) 
    qpAdm_Fixed_sources_summary[nrow(qpAdm_Fixed_sources_summary)+1, ] <- c(as.character(qpAdm_result$weights$target[1]),
                                                                            as.character(qpAdm_result$weights$left[1]),
                                                                            as.character(qpAdm_result$weights$left[2]),
                                                                            as.character(qpAdm_result$weights$left[3]),
                                                                            as.numeric(qpAdm_result$weights$weight[1]),
                                                                            as.numeric(qpAdm_result$weights$weight[2]),
                                                                            as.numeric(qpAdm_result$weights$weight[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$p[1]),
                                                                            as.numeric(qpAdm_result$weights$se[1]),
                                                                            as.numeric(qpAdm_result$weights$se[2]),
                                                                            as.numeric(qpAdm_result$weights$se[3]),
                                                                            as.numeric(qpAdm_result$rankdrop$f4rank[1]),
                                                                            as.numeric(rep))
  }
  
  #' Return qpAdm summaries  
  #return(qpAdm_Fixed_sources_summary)
  
}

#' Return qpAdm summaries  
#return(qpAdm_Fixed_sources_summary)

# Save output
saveRDS(qpAdm_Fixed_sources_summary, file=opt$out)


                                           
