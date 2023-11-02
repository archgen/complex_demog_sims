#' Creating a Hapmap recombination map for msprime RateMap function
library(data.table)
library(ggplot2)

setwd("~/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/")

#' msprime read_hapmap expects a white-space-delimited file with a single header line
#' Columns are: 
#' Chromosome: Can be chr or integer chr10 or 10
#' Position(bp): a physical position (in base pairs)
#' Rate(cM/Mb): recombination rate (in centiMorgans per megabase). 
#'       This gives the constant recom. rate between the physical position in that line (inclusive) and the physical position on the next line (exclusive)
#' Map(cM): a genetic map position (in centiMorgans) 
# ratemap recombination rate == genetic distance span (cM) / physical distance span (converted bp to Mb)

#' NOTE ON RECOMBINATION
#' From msprime manual, "the recombination rate in the base pair segments separating chromosomes to log(2).
#' 
#' The recombination rate is centiMorgans (cMs) per megabase. 
#' 1 cM ~ 1% probability of recombination so for a ~50% probability of recombination we are after 50 cM. 
#' Because we are after a log(2) recombination rate for msprime, we'll thus simulate
#' a log(2)*100 cM recombination probability. 
#' Note that msprime expects the Rate parameter to be cM/Mb (Mb == 1e6). So our rate separating 
#' two chromosomes will be (log(2)*100)*1e6   

#' Source of hapmap genome: https://github.com/adimitromanolakis/geneticMap-GRCh37
hapmap = rbindlist(lapply(1:6, function(x) fread(paste0("./HapMapData/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr", x, ".txt"))))
# Columns: Chromosome Position(bp) Rate(cM/Mb)   Map(cM)
colnames(hapmap) = c("chr", "pos", "rrate", "map")
hapmap[, chr_num := as.numeric(gsub("chr(.+)", "\\1", chr))]

hapmap[, min_pos := min(pos), chr]
hapmap[, min_map := min(map), chr]

hapmap[, new_pos := pos - min_pos + 1]
hapmap[, new_map := map - min_map + (log(2)*100)]


hapmap[, new_rrate := rrate]
hapmap[, new_chr := "chr1"]

hapmap[chr_num == 1, new_pos_cum := new_pos]
hapmap[chr_num == 1, new_map_cum := new_map]

for (i in 2:6) {
  maxPrevPosValue = hapmap[chr_num == i-1, max(new_pos_cum)]
  maxPrevMapValue = hapmap[chr_num == i-1, max(new_map_cum)]
  
  hapmap[chr_num == i, new_pos_cum := maxPrevPosValue + new_pos]
  hapmap[chr_num == i, new_map_cum := maxPrevMapValue + new_map]
  
  #hapmap[chr_num == i-1 & pos == max(pos[chr_num == i-1]), new_rrate := 50*1e6] 
  hapmap[chr_num == i-1 & pos == max(pos[chr_num == i-1]), new_rrate := (log(2)*100)*1e6]
}


filterID = c(TRUE, rep(FALSE, 100))  # To reduced size of plots...
hapmap[, plotID := rep(filterID, length = .N)]  # To reduced size of plots...
hapmap[new_rrate == (log(2)*100)*1e6, plotID := TRUE]

par(mfrow = c(1,2))
hapmap[plotID == TRUE, plot(new_pos_cum, new_map_cum)]
hapmap[plotID == TRUE, plot(new_pos_cum, new_rrate)]

## 
hapmap_out = hapmap[, .(chr_num, new_pos_cum, format(new_rrate, scientific = F), new_map_cum)]
hapmap_out_chr = hapmap[, .(chr, new_pos_cum, format(new_rrate, scientific = F), new_map_cum)]

colnames(hapmap_out) = c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
colnames(hapmap_out_chr) = c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
hapmap_out_chr$Chromosome <- as.character(hapmap_out_chr$Chromosome)

options(scipen = 999)
write.table(hapmap_out, file = "./snakemake/project__complex_msp_simulation/input/rate_map_genetic_map_HapMapII_GRCh37_Chrs1_6__log2recom.txt",
            quote = FALSE, row.names = FALSE)
