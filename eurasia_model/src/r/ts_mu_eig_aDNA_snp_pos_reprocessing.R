#' Re-code the .snp "aDNA" eigenstrat file
#' 
#' Matthew Williams: 
#'  
#' Load in libraries
#.libPaths("/storage/home/mkw5910/.conda/envs/msprime-env/lib/R/library")
library(tidyverse)
library(optparse)
library(data.table)
# SYSTEM ARGUMENTS
option_list = list(
  make_option(c("--sim_snp"),
              type = "character",
              default = NULL,
              help = "geno simulated eigenstrat file",
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
out = opt$out
snp = fread(file = sim_snp, header = F)


#' Chromosome Lengths
chr1_L = 248956422
chr2_L = 242193529
chr3_L = 198295559
chr4_L = 190214555
chr5_L = 181538259
chr6_L = 170805979
chr7_L = 159345973
chr8_L = 145138636
chr9_L = 138394717
chr10_L = 133797422
chr11_L = 135086622
chr12_L = 133275309
chr13_L = 114364328
chr14_L = 107043718
chr15_L = 101991189
chr16_L = 90338345
chr17_L = 83257441
chr18_L = 80373285
chr19_L = 58617616
chr20_L = 64444167
chr21_L = 46709983
chr22_L = 50818468

# Define the chromosome lengths in a list
chromosome_lengths <- c(
  chr1_L, chr2_L, chr3_L, chr4_L, chr5_L,
  chr6_L, chr7_L, chr8_L, chr9_L, chr10_L,
  chr11_L, chr12_L, chr13_L, chr14_L, chr15_L,
  chr16_L, chr17_L, chr18_L, chr19_L, chr20_L,
  chr21_L, chr22_L
)

# Calculate the cumulative sum to get the maximum positions for each chromosome
chr_max <- cumsum(chromosome_lengths)


# Set the directory path
directory_path <- "/Users/mkw5910/Documents/PSU_sims/snakemake/src_PSU_simulations_Projects/complex_ane_demography/roar_collab_out/ts_mu_eig_aDNA_1240K_public_v50/"
# Get the list of files ending with .snp extension
file_list <- list.files(directory_path, pattern = "\\.snp$", full.names = TRUE)

# Print the list of files
print(file_list)
for(f in file_list){
  
  snp = fread(file = f, header = F)
  # Recode the Chromosome position based on the simulated Chromosome lengths
  snp[, V2 := findInterval(V4, c(0,chr_max))]

  #' Recode the phyical position as the current pos - the cumulitative max pos of the preceeding chr
  snp[V2 %between% c(2, 22), new_pos := V4 - chr_max[V2 - 1], by = V2]
  snp[V2 == 1,new_pos := V4]
  
  snp <- snp[,.(V1, V2, V3, new_pos, V5, V6)]
  # Save .snp to file
  fwrite(snp, file = f, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  #' Get a summary of the SNPs counts for each chromosome
  snps_summary = data.table(snp %>% group_by(V2) %>% summarize(counts = n()))
  colnames(snps_summary) <- c("Chr", "Counts")
  # Write the summary to a text file
  write.table(snps_summary, file = paste0(f, "__Chr_snp_counts_summary.txt"), sep = "\t", row.names = FALSE)
}
