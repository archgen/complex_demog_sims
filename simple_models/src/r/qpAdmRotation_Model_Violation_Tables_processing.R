#' Admixtools2 processing qpAdm results & Generating Figures
#' 
#' Matthew Williams
library(data.table);library(dplyr)

## Testing
#path = "/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/project__complex_msp_simulation/"
#data = "out/Simple_Tree_Chr22_rotation/admixtools2/qpAdm_branch/qpAdm_rotation_ts_branch_F2__Simple_Tree_Chr22_rotation_rep_1.rds"
#fileNames = paste0(path, data)

qpAdm_rotate_processing = function(DIR, sim_model, dataType) {
  weightsTable = data.table()
  rankdropTable = data.table()
  popdropTable = data.table()
  sim_Name = sim_model
  dataType=dataType
  #file_DIR = DIR
  fileNames = DIR
  
  #fileNames = dir(file_DIR, pattern = "*rep", full.names = TRUE)
  qpadm_RDS <- lapply(fileNames, readRDS)
  
  for (i in 1:length(qpadm_RDS)) {
    qpAdm_run = qpadm_RDS[[i]]
    weightsTableTemp = rbindlist(lapply(qpAdm_run, function(y) y$weights), idcol = TRUE)
    weightsTableTemp[, admix_model_ID := .id]
    rankdropTableTemp = rbindlist(lapply(qpAdm_run, function(y) y$rankdrop), idcol = TRUE)
    rankdropTableTemp[, admix_model_ID := .id]
    #popdropTableTemp = rbindlist(lapply(qpAdm_run, function(y) y$popdrop), idcol = TRUE, fill=TRUE)
    #popdropTableTemp[, admix_model_ID := .id]
    
    weightsTable = rbind(weightsTable, cbind(sim_model = paste0(sim_Name, "_", dataType), sim = i, weightsTableTemp))
    rankdropTable = rbind(rankdropTable, cbind(sim_model = paste0(sim_Name, "_", dataType), sim = i, rankdropTableTemp))
    #popdropTable = rbind(popdropTable, cbind(sim_model = paste0(sim_Name, "_", dataType), sim = i, popdropTableTemp))
  }
  
  weightsTable[, numLeftPops := .N, .(sim, admix_model_ID, sim_model)]
  rankdropTable[, numLeftPops := .N, .(sim, admix_model_ID, sim_model)]
  #popdropTable[, numLeftPops := .N, .(rep, admix_model_ID, sim_model)]
  # Merge weights and rank tables
  simSummaryTable_tmp1 <- cbind(rankdropTable, weightsTable)
  # subset for columns of interest (sim model, sim , admix_model, numLeftPops, pvalue, left, weight, se, z, chisqdiff, p_nested)
  simSummaryTable_tmp2 <- data.table(simSummaryTable_tmp1[, c(1,2,4,7,9,10,11,12,16,17,18,19,20)])
  ## Generate merged separate datatables for each set of models with c(1,2,3) source populations
  
  ### Creating Summary qpAdm Table when there are models with 3 source pops
  if(nrow(simSummaryTable_tmp2[numLeftPops == 3]) != 0){
    simSummaryTable_S3_tmp <- simSummaryTable_tmp2[numLeftPops == 3]
    # aggregate the samples from the same admixture model into a single line
    simSummaryTable_S3_tmp_tmp2 <- data.table(aggregate(data=simSummaryTable_S3_tmp, cbind(left,weight,p,se)~sim_model+sim+admix_model_ID+numLeftPops+target, FUN=paste))
    simSummaryTable_S3 <- simSummaryTable_S3_tmp_tmp2 %>%
      dplyr::select(sim_model, sim , admix_model_ID, left.V1, left.V2, left.V3, weight.V1, weight.V2, weight.V3,
                    se.V1, se.V2, se.V3,
                    p.V1,
                    target,
                    numLeftPops)
    ### 2 source pops
    simSummaryTable_S2_tmp <- simSummaryTable_tmp2[numLeftPops == 2]
    # aggregate the samples from the same admixture model into a single line
    simSummaryTable_S2_tmp2 <- data.table(aggregate(data=simSummaryTable_S2_tmp, cbind(left,weight,p,se)~sim_model+sim+admix_model_ID+numLeftPops+target, FUN=paste))
    simSummaryTable_S2 <- simSummaryTable_S2_tmp2[, ':=' (left.V3=NA,weight.V3=NA,se.V3=NA)] %>%
      dplyr::select(sim_model, sim , admix_model_ID, left.V1, left.V2, left.V3, weight.V1, weight.V2, weight.V3, se.V1, se.V2, se.V3, p.V1, target, numLeftPops)
    ### 1 source pop
    simSummaryTable_S1_tmp <- simSummaryTable_tmp2[numLeftPops == 1]
    # aggregate the samples from the same admixture model into a single line
    simSummaryTable_S1_tmp2 <- data.table(aggregate(data=simSummaryTable_S1_tmp, cbind(left,weight,p,se)~sim_model+sim+admix_model_ID+numLeftPops+target, FUN=paste))
    setnames(simSummaryTable_S1_tmp2, c("sim_model", "sim", "admix_model_ID", "numLeftPops","target" ,"left.V1", "weight.V1", "p.V1", "se.V1"))
    # Add the columns for models 2 and 3 into the models 1 table
    simSummaryTable_S1 <- simSummaryTable_S1_tmp2[, ':=' (left.V2=NA,left.V3=NA,weight.V2=NA,weight.V3=NA,se.V2=NA,se.V3=NA)] %>%
      dplyr::select(sim_model, sim , admix_model_ID, left.V1, left.V2, left.V3, weight.V1, weight.V2, weight.V3, se.V1, se.V2, se.V3, p.V1, target, numLeftPops)
    
    # Merge the three data.tables into single data table summary of admix models of sources 1-3 & remove old single files
    # Merge S2 and S3 with S1
    simSummaryTable_S1_3 <- rbind(simSummaryTable_S1, simSummaryTable_S2, simSummaryTable_S3)
    # Clean up individual tables
    rm(simSummaryTable_S1)
    rm(simSummaryTable_S2)
    rm(simSummaryTable_S3)
    
    # relabel the column names and set the numerical values to numerical
    setnames(simSummaryTable_S1_3, c("sim_model", "sim", "admix_model_ID", "Source1", "Source2", "Source3", "weight1", "weight2", "weight3","se.weight1", "se.weight2", "se.weight3" ,"pvalue", "target" , "numLeftPops"))
    simSummaryTable_S1_3$pvalue <- as.numeric(simSummaryTable_S1_3$pvalue)
    simSummaryTable_S1_3$weight1 <- as.numeric(simSummaryTable_S1_3$weight1)
    simSummaryTable_S1_3$weight2 <- as.numeric(simSummaryTable_S1_3$weight2)
    simSummaryTable_S1_3$weight3 <- as.numeric(simSummaryTable_S1_3$weight3)
    simSummaryTable_S1_3$se.weight1 <- as.numeric(simSummaryTable_S1_3$se.weight1)
    simSummaryTable_S1_3$se.weight2 <- as.numeric(simSummaryTable_S1_3$se.weight2)
    simSummaryTable_S1_3$se.weight3 <- as.numeric(simSummaryTable_S1_3$se.weight3)
    
    # Define Plausible admixture models
    simSummaryTable_S1_3[pvalue >= 0.01 & numLeftPops == 1, "plausible" := 1]
    simSummaryTable_S1_3[numLeftPops == 2 & pvalue >= 0.01 & weight1 < 1 & weight1 > 0 & weight2 > 0 & weight2 < 1, plausible := 1]
    simSummaryTable_S1_3[numLeftPops ==3 & pvalue >= 0.01 & weight1 < 1 & weight1 > 0 & weight2 > 0 & weight2 < 1 & weight3 > 0 & weight3 < 1, plausible := 1]
    # replace <NA> with 0
    simSummaryTable_S1_3[is.na(simSummaryTable_S1_3$plausible)]$plausible <- 0
    
    ### clean-up temp files
    rm(list = ls(pattern= "tmp"))
    rm(qpadm_RDS)
    
    return(simSummaryTable_S1_3)
    
  } else {
    
    ### 2 source pops
    simSummaryTable_S2_tmp <- simSummaryTable_tmp2[numLeftPops == 2]
    # aggregate the samples from the same admixture model into a single line
    simSummaryTable_S2_tmp_tmp2 <- data.table(aggregate(data=simSummaryTable_S2_tmp, cbind(left,weight,p,se)~sim_model+sim+admix_model_ID+numLeftPops+target, FUN=paste))
    simSummaryTable_S2 <- simSummaryTable_S2_tmp_tmp2 %>%
      dplyr::select(sim_model, sim , admix_model_ID, left.V1, left.V2, weight.V1, weight.V2, 
                    se.V1, se.V2,
                    p.V1,
                    target,
                    numLeftPops)
    ### 1 source pop
    simSummaryTable_S1_tmp <- simSummaryTable_tmp2[numLeftPops == 1]
    # aggregate the samples from the same admixture model into a single line
    simSummaryTable_S1_tmp2 <- data.table(aggregate(data=simSummaryTable_S1_tmp, cbind(left,weight,p,se)~sim_model+sim+admix_model_ID+numLeftPops+target, FUN=paste))
    setnames(simSummaryTable_S1_tmp2, c("sim_model", "sim", "admix_model_ID", "numLeftPops","target" ,"left.V1", "weight.V1", "p.V1", "se.V1"))
    # Add the columns for models 2 and 3 into the models 1 table
    simSummaryTable_S1 <- simSummaryTable_S1_tmp2[, ':=' (left.V2=NA,weight.V2=NA,se.V2=NA)] %>%
      dplyr::select(sim_model, sim , admix_model_ID, left.V1, left.V2, weight.V1, weight.V2, se.V1, se.V2, p.V1, target, numLeftPops)
    
    # Merge the three data.tables into single data table summary of admix models of sources 1-3 & remove old single files
    # Merge S2 and S3 with S1
    simSummaryTable_S1_2 <- rbind(simSummaryTable_S1, simSummaryTable_S2)
    
    # relabel the column names and set the numerical values to numerical
    setnames(simSummaryTable_S1_2, c("sim_model", "sim", "admix_model_ID", "Source1", "Source2", "weight1", "weight2", "se.weight1", "se.weight2","pvalue", "target" , "numLeftPops"))
    simSummaryTable_S1_2$pvalue <- as.numeric(simSummaryTable_S1_2$pvalue)
    simSummaryTable_S1_2$weight1 <- as.numeric(simSummaryTable_S1_2$weight1)
    simSummaryTable_S1_2$weight2 <- as.numeric(simSummaryTable_S1_2$weight2)
    simSummaryTable_S1_2$se.weight1 <- as.numeric(simSummaryTable_S1_2$se.weight1)
    simSummaryTable_S1_2$se.weight2 <- as.numeric(simSummaryTable_S1_2$se.weight2)
    
    # Define Plausible admixture models
    simSummaryTable_S1_2[numLeftPops == 1 & pvalue >= 0.01 & weight1 <= 1 & weight1 >= 0,  "plausible" := 1]
    simSummaryTable_S1_2[numLeftPops == 2 & pvalue >= 0.01 & weight1 <= 1 & weight1 >= 0 & weight2 >= 0 & weight2 <= 1, plausible := 1]
    # replace <NA> with 0
    simSummaryTable_S1_2[is.na(simSummaryTable_S1_2$plausible)]$plausible <- 0
    
    ### clean-up temp files
    rm(list = ls(pattern= "tmp"))
    rm(qpadm_RDS)
    
    return(simSummaryTable_S1_2)
  }
}
