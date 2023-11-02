#' Reading in the complex simulation qpAdm Rotation aDNA and branch
#' files, merging them and generating a summary analysis file with 
#' plausibility criteria evaluated

#' Load libraries.
library(data.table);library(tidyverse)

#' set working directory to output top directory
setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/out/complexDemography_model_A_HapMapII_GRCh37__Chrs21_22/")

#######################
#' qpAdm Rotation Files
branch_qpAdmRotate = fread(file = "./admixtools2/ts_branch_SummaryOutput_qpAdmRotation/ts_branch_qpAdmRotation_SummaryTable__complexDemography_model_A_HapMapII_GRCh37__Chrs21_22.txt")
aDNA_qpAdmRotate = fread(file = "./admixtools2/ts_mu_eig_aDNA_SummaryOutput_qpAdmRotation/ts_mu_eig_aDNA_qpAdmRotation_SummaryTable__complexDemography_model_A_HapMapII_GRCh37__Chrs21_22.txt")

###################
#' Define Variables 
simModel= "complexDemography_model_A_HapMapII_GRCh37__Chrs21_22"
taret_pop = "sLev_IA1"
truemodel_id = 33
#' Visual check
unique(aDNA_qpAdmRotate[admix_model_ID == truemodel_id, c('Source1', 'Source2')]) 
unique(branch_qpAdmRotate[admix_model_ID == truemodel_id, c('Source1', 'Source2')])
Nsims = max(branch_qpAdmRotate$sim)
Nmodels = max(branch_qpAdmRotate$admix_model_ID)

# set the directory path
out_Fig <- "./processed_results/qpAdmRotation/Figures/"

# check if directory exists and create if it doesn't
if (!file.exists(out_Fig)) {
  dir.create(out_Fig, recursive = TRUE)
}

# set the directory path
out_Tab <- "./processed_results/qpAdmRotation/Tables/"

# check if directory exists and create if it doesn't
if (!file.exists(out_Tab)) {
  dir.create(out_Tab, recursive = TRUE)
}


#' Set order of the data tables to match
branch_qpAdmRotate <- branch_qpAdmRotate[order(sim, admix_model_ID)]
aDNA_qpAdmRotate <- aDNA_qpAdmRotate[order(sim, admix_model_ID)]
identical(branch_qpAdmRotate[, sim, admix_model_ID], aDNA_qpAdmRotate[, sim, admix_model_ID])
aDNA_qpAdmRotate[, Row_Number := .I]
branch_qpAdmRotate[, Row_Number := .I]

#' Check that Single-Source models have the same Soruce1 pop
identical(branch_qpAdmRotate[numLeftPops == 1, Source1], aDNA_qpAdmRotate[numLeftPops == 1, Source1])

#' Get the Source1, Source2 and Source3 populations 
#' to match between the branch and aDNA datasets

#' For 2-source models 
branch_qpAdmRotate_2Sources = branch_qpAdmRotate[numLeftPops == 2]
aDNA_qpAdmRotate_2Sources = aDNA_qpAdmRotate[numLeftPops == 2]
for(i in 1:nrow(branch_qpAdmRotate_2Sources)){
  #' Of the three source pops in the aDNA DT, which position does that branch DT source pop match ..? .. Save that position
  branch_S1pop_aDNApos = which(aDNA_qpAdmRotate_2Sources[i, c('Source1', 'Source2')] %in% branch_qpAdmRotate_2Sources[i, 'Source1'])
  branch_S2pop_aDNApos = which(aDNA_qpAdmRotate_2Sources[i, c('Source1', 'Source2')] %in% branch_qpAdmRotate_2Sources[i, 'Source2'])

  #' Create a new position for each respective branch DT source (1,2 and 3) that was obtained from the matching above. 
  s1_new_pos = paste0("Source", as.character(branch_S1pop_aDNApos))
  s2_new_pos = paste0("Source", as.character(branch_S2pop_aDNApos))

  #' Paste the original branch DT source pop into it's new position which will match the position that it appears in in the aDNA DT 
  branch_qpAdmRotate_2Sources[i, paste0(s1_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('Source1')]$Source1]
  branch_qpAdmRotate_2Sources[i, paste0(s2_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('Source2')]$Source2]

  #' do the same as above for the weights
  w1_new_pos = paste0("weight", as.character(branch_S1pop_aDNApos))
  w2_new_pos = paste0("weight", as.character(branch_S2pop_aDNApos))

  branch_qpAdmRotate_2Sources[i, paste0(w1_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('weight1')]$weight1]
  branch_qpAdmRotate_2Sources[i, paste0(w2_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('weight2')]$weight2]

  #' do the same as above for the se
  se1_new_pos = paste0("se.", as.character(branch_S1pop_aDNApos))
  se2_new_pos = paste0("se.", as.character(branch_S2pop_aDNApos))

  branch_qpAdmRotate_2Sources[i, paste0(se1_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('se.weight1')]$se.weight1]
  branch_qpAdmRotate_2Sources[i, paste0(se2_new_pos, "_aDNApos") := branch_qpAdmRotate_2Sources[i, c('se.weight2')]$se.weight2]

}
#' Visual check
identical(branch_qpAdmRotate_2Sources[, Source2_aDNApos], aDNA_qpAdmRotate_2Sources[, Source2])
identical(branch_qpAdmRotate_2Sources[, Source1_aDNApos], aDNA_qpAdmRotate_2Sources[, Source1])
branch_qpAdmRotate_2Sources[2790] %>% select(Source1_aDNApos,Source2_aDNApos)
aDNA_qpAdmRotate[numLeftPops == 2][2790] %>% select(Source1, Source2)

#' Re-fil the original dataset with the updated re-sorted source positions and values
branch_qpAdmRotate[numLeftPops == 2, c('Source1', 'Source2',  
                                       'weight1', 'weight2', 
                                       'se.weight1', 'se.weight2')] <- branch_qpAdmRotate_2Sources[, c('Source1_aDNApos', 'Source2_aDNApos',  
                                                                                                                     'weight1_aDNApos', 'weight2_aDNApos',
                                                                                                                     'se.1_aDNApos', 'se.2_aDNApos')]

identical(branch_qpAdmRotate[numLeftPops == 2, Source1], aDNA_qpAdmRotate[numLeftPops == 2, Source1])


#' For 3-source models 
branch_qpAdmRotate_3Sources = branch_qpAdmRotate[numLeftPops == 3]
aDNA_qpAdmRotate_3Sources = aDNA_qpAdmRotate[numLeftPops == 3]
for(i in 1:nrow(branch_qpAdmRotate_3Sources)){
  #' Of the three source pops in the aDNA DT, which position does that branch DT source pop match ..? .. Save that position
  branch_S1pop_aDNApos = which(aDNA_qpAdmRotate_3Sources[i, c('Source1', 'Source2', 'Source3')] %in% branch_qpAdmRotate_3Sources[i, 'Source1'])
  branch_S2pop_aDNApos = which(aDNA_qpAdmRotate_3Sources[i, c('Source1', 'Source2', 'Source3')] %in% branch_qpAdmRotate_3Sources[i, 'Source2'])
  branch_S3pop_aDNApos = which(aDNA_qpAdmRotate_3Sources[i, c('Source1', 'Source2', 'Source3')] %in% branch_qpAdmRotate_3Sources[i, 'Source3'])
  
  #' Create a new position for each respective branch DT source (1,2 and 3) that was obtained from the matching above. 
  s1_new_pos = paste0("Source", as.character(branch_S1pop_aDNApos))
  s2_new_pos = paste0("Source", as.character(branch_S2pop_aDNApos))
  s3_new_pos = paste0("Source", as.character(branch_S3pop_aDNApos))
  
  #' Paste the original branch DT source pop into it's new position which will match the position that it appears in in the aDNA DT 
  branch_qpAdmRotate_3Sources[i, paste0(s1_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('Source1')]$Source1]
  branch_qpAdmRotate_3Sources[i, paste0(s2_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('Source2')]$Source2]
  branch_qpAdmRotate_3Sources[i, paste0(s3_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('Source3')]$Source3]

  
  #' do the same as above for the weights
  w1_new_pos = paste0("weight", as.character(branch_S1pop_aDNApos))
  w2_new_pos = paste0("weight", as.character(branch_S2pop_aDNApos))
  w3_new_pos = paste0("weight", as.character(branch_S3pop_aDNApos))
  
  branch_qpAdmRotate_3Sources[i, paste0(w1_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('weight1')]$weight1]
  branch_qpAdmRotate_3Sources[i, paste0(w2_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('weight2')]$weight2]
  branch_qpAdmRotate_3Sources[i, paste0(w3_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('weight3')]$weight3]
  
  #' do the same as above for the se
  se1_new_pos = paste0("se.", as.character(branch_S1pop_aDNApos))
  se2_new_pos = paste0("se.", as.character(branch_S2pop_aDNApos))
  se3_new_pos = paste0("se.", as.character(branch_S3pop_aDNApos))
  
  branch_qpAdmRotate_3Sources[i, paste0(se1_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('se.weight1')]$se.weight1]
  branch_qpAdmRotate_3Sources[i, paste0(se2_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('se.weight2')]$se.weight2]
  branch_qpAdmRotate_3Sources[i, paste0(se3_new_pos, "_aDNApos") := branch_qpAdmRotate_3Sources[i, c('se.weight3')]$se.weight3]
  
}
#' Visual Check
branch_qpAdmRotate_3Sources[50] %>% select(Source1, weight1, se.weight1, Source2, weight2, se.weight2, Source3, weight3, se.weight3)
aDNA_qpAdmRotate[numLeftPops == 3][50] %>% select(Source1, weight1, Source2,  weight2, Source3, weight3)
branch_qpAdmRotate_3Sources[50] %>% select(Source1_aDNApos, weight1_aDNApos, se.1_aDNApos, 
                                           Source2_aDNApos, weight2_aDNApos,  se.2_aDNApos, 
                                           Source3_aDNApos, weight3_aDNApos,  se.3_aDNApos)

#' Re-fil the original dataset with the updated re-sorted source positions and values
branch_qpAdmRotate[numLeftPops == 3, c('Source1', 'Source2', 'Source3', 
                                       'weight1', 'weight2', 'weight3', 
                                       'se.weight1', 'se.weight2', 'se.weight3')] <- branch_qpAdmRotate_3Sources[, c('Source1_aDNApos', 'Source2_aDNApos', 'Source3_aDNApos', 
                                                                                                                     'weight1_aDNApos', 'weight2_aDNApos', 'weight3_aDNApos', 
                                                                                                                     'se.1_aDNApos', 'se.2_aDNApos', 'se.3_aDNApos')]

identical(branch_qpAdmRotate[numLeftPops == 3, Source1], aDNA_qpAdmRotate[numLeftPops == 3, Source1])
identical(branch_qpAdmRotate[numLeftPops == 3, Source2], aDNA_qpAdmRotate[numLeftPops == 3, Source2])
identical(branch_qpAdmRotate[numLeftPops == 3, Source3], aDNA_qpAdmRotate[numLeftPops == 3, Source3])

#' All source populations equal? 
identical(branch_qpAdmRotate[, Source1], aDNA_qpAdmRotate[, Source1])
identical(branch_qpAdmRotate[, Source2], aDNA_qpAdmRotate[, Source2])
identical(branch_qpAdmRotate[, Source3], aDNA_qpAdmRotate[, Source3])



### ********************************************** ###
### ***** PLAUSIBILITY CRITERIA & EVALUATION ***** ###


#' Generate Plausibility Conditions
for(i in seq(1:nrow(aDNA_qpAdmRotate))) {
  if(aDNA_qpAdmRotate[i, pvalue >= 0.05 & plausible == 1]){
    aDNA_qpAdmRotate[i, "plausible_05" := 1]
  } else{
    aDNA_qpAdmRotate[i, "plausible_05" := 0]
  }
}
for(i in seq(1:nrow(branch_qpAdmRotate))) {
  if(branch_qpAdmRotate[i, pvalue >= 0.05 & plausible == 1]){
    branch_qpAdmRotate[i, "plausible_05" := 1]
  } else{
    branch_qpAdmRotate[i, "plausible_05" := 0]
  }
}


#' QTP p-value > 0.05 & weights [0:1]
#' Branch Lengths
for (i in 1:Nsims){
  if(branch_qpAdmRotate[sim == i & admix_model_ID == truemodel_id]$plausible_05 != 1){ 
    # If the true model is not plausible at 0.05 criteria: 
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ) minus 1). Max value = 0 EQ{(1 - (0 / 20)) - 1}, Min value = -1 EQ{(1 - (20 / 20)) - 1} : N = 20. 
    branch_qpAdmRotate[sim == i, 
                               'QTP_p005_weight01':= (1 - ((nrow(branch_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible_05 == 1]) / (Nmodels - 1)))) - 1]
    # If the true model is plausible at 0.05 criteria:
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ). Max value = 1 EQ{1 - (0 / 20)}, Min value = 0 EQ{1 - (20 / 20)} : N = 20. 
  } else {
    branch_qpAdmRotate[sim == i, 
                               'QTP_p005_weight01':= 1 - ((nrow(branch_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible_05 == 1]) / (Nmodels - 1)))]
  }
}
#' aDNA sim 
for (i in 1:Nsims){
  if(aDNA_qpAdmRotate[sim == i & admix_model_ID == truemodel_id]$plausible_05 != 1){ 
    # If the true model is not plausible at 0.05 criteria: 
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ) minus 1). Max value = 0 EQ{(1 - (0 / 20)) - 1}, Min value = -1 EQ{(1 - (20 / 20)) - 1} : N = 20. 
    aDNA_qpAdmRotate[sim == i, 
                             'QTP_p005_weight01':= (1 - ((nrow(aDNA_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible_05 == 1]) / (Nmodels - 1)))) - 1]
    # If the true model is plausible at 0.05 criteria:
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ). Max value = 1 EQ{1 - (0 / 20)}, Min value = 0 EQ{1 - (20 / 20)} : N = 20. 
  } else {
    aDNA_qpAdmRotate[sim == i, 
                             'QTP_p005_weight01':= 1 - ((nrow(aDNA_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible_05 == 1]) / (Nmodels - 1)))]
  }
}

#' Binary outcome for QTP WITH Stringent p-value & relaxed weights
aDNA_qpAdmRotate$QTP_p005_weight01_binary = ifelse(aDNA_qpAdmRotate$QTP_p005_weight01 != 1, 0, 1)
branch_qpAdmRotate$QTP_p005_weight01_binary = ifelse(branch_qpAdmRotate$QTP_p005_weight01 != 1, 0, 1)


#' QTP p-value > 0.01 & weights [0:1]
#' Branch Lengths
for (i in 1:Nsims){
  if(branch_qpAdmRotate[sim == i & admix_model_ID == truemodel_id]$plausible != 1){ 
    # If the true model is not plausible at 0.05 criteria: 
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ) minus 1). Max value = 0 EQ{(1 - (0 / 20)) - 1}, Min value = -1 EQ{(1 - (20 / 20)) - 1} : N = 20. 
    branch_qpAdmRotate[sim == i, 
                               'QTP_p001_weight01':= (1 - ((nrow(branch_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible == 1]) / (Nmodels - 1)))) - 1]
    # If the true model is plausible at 0.05 criteria:
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ). Max value = 1 EQ{1 - (0 / 20)}, Min value = 0 EQ{1 - (20 / 20)} : N = 20. 
  } else {
    branch_qpAdmRotate[sim == i, 
                               'QTP_p001_weight01':= 1 - ((nrow(branch_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible == 1]) / (Nmodels - 1)))]
  }
}
#' aDNA sim 
for (i in 1:Nsims){
  if(aDNA_qpAdmRotate[sim == i & admix_model_ID == truemodel_id]$plausible != 1){ 
    # If the true model is not plausible at 0.05 criteria: 
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ) minus 1). Max value = 0 EQ{(1 - (0 / 20)) - 1}, Min value = -1 EQ{(1 - (20 / 20)) - 1} : N = 20. 
    aDNA_qpAdmRotate[sim == i, 
                             'QTP_p001_weight01':= (1 - ((nrow(aDNA_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible == 1]) / (Nmodels - 1)))) - 1]
    # If the true model is plausible at 0.05 criteria:
    # QTP = ((One minus the number of false models qpAdm deems plausible / 1 - Total number of models ). Max value = 1 EQ{1 - (0 / 20)}, Min value = 0 EQ{1 - (20 / 20)} : N = 20. 
  } else {
    aDNA_qpAdmRotate[sim == i, 
                             'QTP_p001_weight01':= 1 - ((nrow(aDNA_qpAdmRotate[sim == i & admix_model_ID != truemodel_id & plausible == 1]) / (Nmodels - 1)))]
  }
}

#' Binary outcome for QTP WITH Stringent p-value & relaxed weights
aDNA_qpAdmRotate$QTP_p001_weight01_binary = ifelse(aDNA_qpAdmRotate$QTP_p001_weight01 != 1, 0, 1)
branch_qpAdmRotate$QTP_p001_weight01_binary = ifelse(branch_qpAdmRotate$QTP_p001_weight01 != 1, 0, 1)

#' Visual Inspection
branch_qpAdmRotate %>% filter(sim == 58) %>% select(admix_model_ID, weight1, weight2,  weight3, 
                                                            pvalue, plausible_05, QTP_p005_weight01, QTP_p005_weight01_binary,
                                                            plausible, QTP_p001_weight01, QTP_p001_weight01_binary)
aDNA_qpAdmRotate %>% filter(sim == 58) %>% select(admix_model_ID, weight1, weight2,  weight3, 
                                                          pvalue, plausible_05, QTP_p005_weight01, QTP_p005_weight01_binary,
                                                          plausible, QTP_p001_weight01, QTP_p001_weight01_binary)


#' Average QTP Binary thresholds 
mean(branch_qpAdmRotate[admix_model_ID == 1]$QTP_p005_weight01)
mean(aDNA_qpAdmRotate[admix_model_ID == 1]$QTP_p005_weight01)

mean(branch_qpAdmRotate[admix_model_ID == 1]$QTP_p001_weight01)
mean(aDNA_qpAdmRotate[admix_model_ID == 1]$QTP_p001_weight01)

mean(branch_qpAdmRotate[admix_model_ID == 1]$QTP_p005_weight01_binary)
mean(aDNA_qpAdmRotate[admix_model_ID == 1]$QTP_p005_weight01_binary)

mean(branch_qpAdmRotate[admix_model_ID == 1]$QTP_p001_weight01_binary)
mean(aDNA_qpAdmRotate[admix_model_ID == 1]$QTP_p001_weight01_binary)




#' Frequency of sims with p > 0.05 and weights outside [0:1]
branch_qpAdmRotate[, fq_p005_weight_outside01 := sum((pvalue >= 0.05 & (weight1 > 1 | weight1 < 0 & weight2 > 1 | weight2 < 0))) / Nreps, by = "admix_model_ID"] 
aDNA_qpAdmRotate[, fq_p005_weight_outside01 := sum((pvalue >= 0.05 & (weight1 > 1 | weight1 < 0 & weight2 > 1 | weight2 < 0))) / Nreps, by = "admix_model_ID"] 

#' Visual Check
nrow(branch_qpAdmRotate[pvalue >= 0.05 & (weight1 > 1 | weight1 < 0 & weight2 > 1 | weight2 < 0) & admix_model_ID == 15]) / Nsims
branch_qpAdmRotate %>% filter(admix_model_ID == 15) %>% select(fq_p005_weight_outside01)
nrow(aDNA_qpAdmRotate[pvalue >= 0.05 & (weight1 > 1 | weight1 < 0 & weight2 > 1 | weight2 < 0) & admix_model_ID == 15]) / Nsims
aDNA_qpAdmRotate %>% filter(admix_model_ID == 15) %>% select(fq_p005_weight_outside01)




#' Tables of Summary Results
#' 
#' aDNA
sourcepops <- unique(aDNA_qpAdmRotate$Source1)
freq_plausible_Sources_Summary_aDNA = data.table("Source" = as.character(),
                                                 "Freq_in_models_p005_weights01" = as.numeric(),
                                                 "Freq_in_models_p001_weights01" = as.numeric(),
                                                 "Freq_in_models_p005_weights2se01" = as.numeric(),
                                                 "Freq_in_models_p001_weights2se01" = as.numeric())
for(i in 1:length(sourcepops)){
  Source <- sourcepops[i]
  freq_plausible_Sources_aDNA <- data.table(Source, "Freq_in_models_p005_weights01" = nrow(aDNA_qpAdmRotate %>% 
                                                                                             filter(plausible == 1 & pvalue >= 0.05 & 
                                                                                                      paste0(Source1, Source2) %like% sourcepops[i])) / max(aDNA_qpAdmRotate$admix_model_ID),
                                            
                                            "Freq_in_models_p001_weights01" = nrow(aDNA_qpAdmRotate %>% 
                                                                                     filter(plausible == 1 & 
                                                                                              paste0(Source1, Source2) %like% sourcepops[i])) / max(aDNA_qpAdmRotate$admix_model_ID),
                                            
                                            "Freq_in_models_p005_weights2se01" = nrow(aDNA_qpAdmRotate %>% 
                                                                                        filter(plausible == 1 & numLeftPops == 2 & pvalue >= 0.05 & (weight1 >= 0 + 2* as.numeric(se.weight1) & weight1 <= 1 - 2* as.numeric(se.weight1)) &
                                                                                                 (weight2 >= 0 + 2* as.numeric(se.weight2) & weight2 <= 1 - 2* as.numeric(se.weight2)) &
                                                                                                 paste0(Source1, Source2) %like% sourcepops[i])) / max(aDNA_qpAdmRotate$admix_model_ID),
                                            
                                            "Freq_in_models_p001_weights2se01" = nrow(aDNA_qpAdmRotate %>% 
                                                                                        filter(plausible == 1 & numLeftPops == 2 & (weight1 >= 0 + 2* as.numeric(se.weight1) & weight1 <= 1 - 2* as.numeric(se.weight1)) &
                                                                                                 (weight2 >= 0 + 2* as.numeric(se.weight2) & weight2 <= 1 - 2* as.numeric(se.weight2)) &
                                                                                                 paste0(Source1, Source2) %like% sourcepops[i])) / max(aDNA_qpAdmRotate$admix_model_ID))
  
  freq_plausible_Sources_Summary_aDNA <- rbind(freq_plausible_Sources_Summary_aDNA, freq_plausible_Sources_aDNA)
  
}

#' Branch
sourcepops <- unique(branch_qpAdmRotate$Source1)
freq_plausible_Sources_Summary_branch = data.table("Source" = as.character(),
                                                 "Freq_in_models_p005_weights01" = as.numeric(),
                                                 "Freq_in_models_p001_weights01" = as.numeric(),
                                                 "Freq_in_models_p005_weights2se01" = as.numeric(),
                                                 "Freq_in_models_p001_weights2se01" = as.numeric())
# for(i in 1:length(sourcepops)){
#   Source <- sourcepops[i]
#   freq_plausible_Sources_branch <- data.table(Source, "Freq_in_models_p005_weights01" = nrow(branch_qpAdmRotate %>% 
#                                                                                              filter(plausible == 1 & pvalue >= 0.05 & 
#                                                                                                       paste0(Source1, Source2) == sourcepops[i])) / max(branch_qpAdmRotate$admix_model_ID),
#                                             
#                                             "Freq_in_models_p001_weights01" = nrow(branch_qpAdmRotate %>% 
#                                                                                      filter(plausible == 1 & 
#                                                                                               paste0(Source1, Source2) == sourcepops[i])) / max(branch_qpAdmRotate$admix_model_ID),
#                                             
#                                             "Freq_in_models_p005_weights2se01" = nrow(branch_qpAdmRotate %>% 
#                                                                                         filter(plausible == 1 & numLeftPops == 2 & pvalue >= 0.05 & (weight1 >= 0 + 2* as.numeric(se.weight1) & weight1 <= 1 - 2* as.numeric(se.weight1)) &
#                                                                                                  (weight2 >= 0 + 2* as.numeric(se.weight2) & weight2 <= 1 - 2* as.numeric(se.weight2)) &
#                                                                                                  paste0(Source1, Source2) == sourcepops[i])) / max(branch_qpAdmRotate$admix_model_ID),
#                                             
#                                             "Freq_in_models_p001_weights2se01" = nrow(branch_qpAdmRotate %>% 
#                                                                                         filter(plausible == 1 & numLeftPops == 2 & (weight1 >= 0 + 2* as.numeric(se.weight1) & weight1 <= 1 - 2* as.numeric(se.weight1)) &
#                                                                                                  (weight2 >= 0 + 2* as.numeric(se.weight2) & weight2 <= 1 - 2* as.numeric(se.weight2)) &
#                                                                                                  paste0(Source1, Source2) == sourcepops[i])) / max(branch_qpAdmRotate$admix_model_ID))
#   
#   freq_plausible_Sources_Summary_branch <- rbind(freq_plausible_Sources_Summary_branch, freq_plausible_Sources_branch)
#   
# }
#' GLANCE
freq_plausible_Sources_Summary_branch
freq_plausible_Sources_Summary_aDNA





###################################################
#' Subset data on 2 Source Models & Generate Labels 
#####################
#' Branch Summary Table 
qpAdm_summary_table_branch_2Source <- branch_qpAdmRotate %>% filter(numLeftPops == 2)
qpAdm_summary_table_branch_2Source$plausible[qpAdm_summary_table_branch_2Source$plausible==1] <- "True"
qpAdm_summary_table_branch_2Source$plausible[qpAdm_summary_table_branch_2Source$plausible==0] <- "False"
qpAdm_summary_table_branch_2Source$plausible_05[qpAdm_summary_table_branch_2Source$plausible_05==1] <- "True"
qpAdm_summary_table_branch_2Source$plausible_05[qpAdm_summary_table_branch_2Source$plausible_05==0] <- "False"
####
qpAdm_summary_table_branch_2Source[, meanPlausible_01 := sum(plausible == "True") / .N, by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, meanPlausible_05 := sum(plausible == "True" & pvalue >= 0.05) / .N, by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Label := as.character(paste0("Model ID=",admix_model_ID)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, PopLabel := as.character(paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Xlab_01 := (paste0("Model ID=",admix_model_ID, "\n", "\U011D", "=", meanPlausible_01)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Xlab_01_FullLabel := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", meanPlausible_01)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Xlab_05 := (paste0("Model ID=",admix_model_ID, "\n", "\U011D", "=", meanPlausible_05)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Xlab_05_FullLabel := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", meanPlausible_05)), by = admix_model_ID]
qpAdm_summary_table_branch_2Source[, Xlab_fq_pWeights := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", fq_p005_weight_outside01)), by = admix_model_ID]
# FACTOR
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, admix_model_ID = factor(admix_model_ID, unique(admix_model_ID)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, Xlab_01 = factor(Xlab_01, unique(Xlab_01)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, Xlab_05 = factor(Xlab_05, unique(Xlab_05)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, PopLabel = factor(PopLabel, unique(PopLabel)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, Xlab_01_FullLabel = factor(Xlab_01_FullLabel, unique(Xlab_01_FullLabel)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, Xlab_05_FullLabel = factor(Xlab_05_FullLabel, unique(Xlab_05_FullLabel)))
qpAdm_summary_table_branch_2Source <- mutate(qpAdm_summary_table_branch_2Source, Xlab_fq_pWeights = factor(Xlab_fq_pWeights, unique(Xlab_fq_pWeights)))


#####################
#####################
#' aDNA Summary Table 
qpAdm_summary_table_aDNA_2Source <- aDNA_qpAdmRotate %>% filter(numLeftPops == 2)
qpAdm_summary_table_aDNA_2Source$plausible[qpAdm_summary_table_aDNA_2Source$plausible==1] <- "True"
qpAdm_summary_table_aDNA_2Source$plausible[qpAdm_summary_table_aDNA_2Source$plausible==0] <- "False"
qpAdm_summary_table_aDNA_2Source$plausible_05[qpAdm_summary_table_aDNA_2Source$plausible_05==1] <- "True"
qpAdm_summary_table_aDNA_2Source$plausible_05[qpAdm_summary_table_aDNA_2Source$plausible_05==0] <- "False"
####
qpAdm_summary_table_aDNA_2Source[, meanPlausible_01 := sum(plausible == "True") / .N, by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, meanPlausible_05 := sum(plausible == "True" & pvalue >= 0.05) / .N, by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Label := as.character(paste0("Model ID=",admix_model_ID)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, PopLabel := as.character(paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Xlab_01 := (paste0("Model ID=",admix_model_ID, "\n", "\U011D", "=", meanPlausible_01)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Xlab_01_FullLabel := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", meanPlausible_01)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Xlab_05 := (paste0("Model ID=",admix_model_ID, "\n", "\U011D", "=", meanPlausible_05)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Xlab_05_FullLabel := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", meanPlausible_05)), by = admix_model_ID]
qpAdm_summary_table_aDNA_2Source[, Xlab_fq_pWeights := (paste0("Model ID=",admix_model_ID, "\n" ,"P1=", Source1, "\n", "P2=", Source2, "\n", "\U011D", "=", fq_p005_weight_outside01)), by = admix_model_ID]
# FACTOR
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, admix_model_ID = factor(admix_model_ID, unique(admix_model_ID)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, Xlab_01 = factor(Xlab_01, unique(Xlab_01)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, Xlab_05 = factor(Xlab_05, unique(Xlab_05)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, PopLabel = factor(PopLabel, unique(PopLabel)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, Xlab_01_FullLabel = factor(Xlab_01_FullLabel, unique(Xlab_01_FullLabel)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, Xlab_05_FullLabel = factor(Xlab_05_FullLabel, unique(Xlab_05_FullLabel)))
qpAdm_summary_table_aDNA_2Source <- mutate(qpAdm_summary_table_aDNA_2Source, Xlab_fq_pWeights = factor(Xlab_fq_pWeights, unique(Xlab_fq_pWeights)))



#' Generate Long-Format Table for Frequency Plausible Figures
#####
# 0.01 Plausible Frequency
qpAdm_summary_table_aDNA_RepSummary_01 <- dcast(aDNA_qpAdmRotate %>% filter(numLeftPops %in% c(1,2)), 
                                                Source1+Source2+admix_model_ID~sim, 
                                                fun=list(sum), 
                                                value.var="plausible")

qpAdm_summary_table_branch_RepSummary_01 <- dcast(branch_qpAdmRotate %>% filter(numLeftPops %in% c(1,2)), 
                                                  Source1+Source2+admix_model_ID~sim, 
                                                  fun=list(sum), 
                                                  value.var="plausible")

# Set the order of the branch to match the aDNA sim, ordered on the admix_model_ID column
setorder(qpAdm_summary_table_branch_RepSummary_01[, .r := order(match(qpAdm_summary_table_aDNA_RepSummary_01$admix_model_ID,qpAdm_summary_table_branch_RepSummary_01$admix_model_ID)
)], .r)[, .r := NULL]

# Then set the Source 1 and Source 2 column names of the branch to be that from the sim aDNA
qpAdm_summary_table_branch_RepSummary_01 = cbind(qpAdm_summary_table_aDNA_RepSummary_01 %>% select(Source1, Source2), qpAdm_summary_table_branch_RepSummary_01[ ,`:=`(Source1 = NULL, Source2 = NULL)] )

# Calculate the mean number of plausible models
qpAdm_summary_table_branch_RepSummary_01[,"mean_plausible_01" := rowMeans(qpAdm_summary_table_branch_RepSummary_01[, 4:ncol(qpAdm_summary_table_branch_RepSummary_01)])]
qpAdm_summary_table_aDNA_RepSummary_01[,"mean_plausible_01" := rowMeans(qpAdm_summary_table_aDNA_RepSummary_01[, 4:ncol(qpAdm_summary_table_aDNA_RepSummary_01)])]

# 0.05 Plausible Frequency
qpAdm_summary_table_aDNA_RepSummary_05 <- dcast(aDNA_qpAdmRotate %>% filter(numLeftPops %in% c(1,2)), 
                                                Source1+Source2+admix_model_ID~sim, 
                                                fun=list(sum), 
                                                value.var="plausible_05")
qpAdm_summary_table_branch_RepSummary_05 <- dcast(branch_qpAdmRotate %>% filter(numLeftPops %in% c(1,2)), 
                                                  Source1+Source2+admix_model_ID~sim, 
                                                  fun=list(sum), 
                                                  value.var="plausible_05")

setorder(qpAdm_summary_table_branch_RepSummary_05[, .r := order(match(qpAdm_summary_table_aDNA_RepSummary_05$admix_model_ID,qpAdm_summary_table_branch_RepSummary_05$admix_model_ID)
)], .r)[, .r := NULL]

qpAdm_summary_table_branch_RepSummary_05 = cbind(qpAdm_summary_table_aDNA_RepSummary_05 %>% select(Source1, Source2), qpAdm_summary_table_branch_RepSummary_05[ ,`:=`(Source1 = NULL, Source2 = NULL)] )


qpAdm_summary_table_branch_RepSummary_05[,"mean_plausible_05" := rowMeans(qpAdm_summary_table_branch_RepSummary_05[, 4:ncol(qpAdm_summary_table_branch_RepSummary_05)])]
qpAdm_summary_table_aDNA_RepSummary_05[,"mean_plausible_05" := rowMeans(qpAdm_summary_table_aDNA_RepSummary_05[, 4:ncol(qpAdm_summary_table_aDNA_RepSummary_05)])]








### ********************************************** ###
### ***** OUTPUT SUMMARY TABLES TO FILE ***** ###

### ***** aDNA Summary Files
saveRDS(qpAdm_summary_table_aDNA_2Source, file = paste0(out_Tab, "qpAdm_summary_table_aDNA_2Source", ".rds"))
saveRDS(aDNA_qpAdmRotate, file = paste0(out_Tab, "aDNA_qpAdmRotate", ".rds"))
saveRDS(qpAdm_summary_table_aDNA_RepSummary_05, file = paste0(out_Tab, "qpAdm_summary_table_aDNA_RepSummary_05", ".rds"))
saveRDS(qpAdm_summary_table_aDNA_RepSummary_01, file = paste0(out_Tab, "qpAdm_summary_table_aDNA_RepSummary_01", ".rds"))




### ***** Branch Summary Files
saveRDS(qpAdm_summary_table_branch_2Source, file = paste0(out_Tab, "qpAdm_summary_table_branch_2Source", ".rds"))
saveRDS(branch_qpAdmRotate, file = paste0(out_Tab, "branch_qpAdmRotate", ".rds"))
saveRDS(qpAdm_summary_table_branch_RepSummary_05, file = paste0(out_Tab, "qpAdm_summary_table_branch_RepSummary_05", ".rds"))
saveRDS(qpAdm_summary_table_branch_RepSummary_01, file = paste0(out_Tab, "qpAdm_summary_table_branch_RepSummary_01", ".rds"))










#' Column bind the two data sets
#' first set the column names to include their data-source label
# Add prefix "branch__" to column names in dt1
# setnames(branch_qpAdmRotate, paste0("branch__", names(branch_qpAdmRotate)))
# # Add prefix "aDNA__" to column names in dt2
# setnames(aDNA_qpAdmRotate, paste0("aDNA__", names(aDNA_qpAdmRotate)))
# # cbind the data.tables
# qpAdmRotate_merged = cbind(branch_qpAdmRotate, aDNA_qpAdmRotate)
