#' Generate qpAdm Rotation Figures for branch and aDNA simulations
#' 
#' 
#' 
#' Load in libraries
library(data.table);library(ggplot2);library(dplyr);library('gtools');library(ggrepel);library(cowplot);library(grid)

###################
#' set working directory to output top directory
setwd("/Users/mkw5910/Library/CloudStorage/Box-Box/PSU/Complex_Admixture_Histories/snakemake/src_PSU_simulations_Projects/complex_demographic_model/out/complexDemography_model_A_HapMapII_GRCh37__Chrs21_22/")

###################
#' Read in Summary qpAdm files
### ***** aDNA Summary Files
qpAdm_summary_table_aDNA_2Source = readRDS(file = paste0(out_Tab, "qpAdm_summary_table_aDNA_2Source", ".rds"))
aDNA_qpAdmRotate = readRDS(file = paste0(out_Tab, "aDNA_qpAdmRotate", ".rds"))
qpAdm_summary_table_aDNA_RepSummary_05 = readRDS(file = paste0(out_Tab, "qpAdm_summary_table_aDNA_RepSummary_05", ".rds"))
qpAdm_summary_table_aDNA_RepSummary_01 = readRDS(file = paste0(out_Tab, "qpAdm_summary_table_aDNA_RepSummary_01", ".rds"))

### ***** Branch Summary Files
qpAdm_summary_table_branch_2Source = readRDS(file = paste0(out_Tab, "qpAdm_summary_table_branch_2Source", ".rds"))
branch_qpAdmRotate = readRDS(file = paste0(out_Tab, "branch_qpAdmRotate", ".rds"))
qpAdm_summary_table_branch_RepSummary_05 = readRDS(file = paste0(out_Tab, "qpAdm_summary_table_branch_RepSummary_05", ".rds"))
qpAdm_summary_table_branch_RepSummary_01= readRDS(file = paste0(out_Tab, "qpAdm_summary_table_branch_RepSummary_01", ".rds"))


###################
#' Define Variables 
simModel= "complexDemography_model_A_HapMapII_GRCh37__Chrs21_22"
taret_pop = "sLev_IA1"
truemodel_id = 33
Nsims = max(branch_qpAdmRotate$sim)
Nmodels = max(branch_qpAdmRotate$admix_model_ID)
out_Fig = "./processed_results/qpAdmRotation/Figures/"
# check if directory exists and create if it doesn't
if (!file.exists(out_Fig)) {
  dir.create(out_Fig, recursive = TRUE)
}

out_Tab = "./processed_results/qpAdmRotation/Tables/"
if (!file.exists(out_Tab)) {
  dir.create(out_Tab, recursive = TRUE)
}



#' Visual check
unique(aDNA_qpAdmRotate[admix_model_ID == truemodel_id, c('Source1', 'Source2')]) 
unique(branch_qpAdmRotate[admix_model_ID == truemodel_id, c('Source1', 'Source2')])




# Create a matrix of data
my_matrix <- matrix(rnorm(100), ncol = 10)

# Convert the matrix to a data frame
df <- reshape2::melt(my_matrix)

# Create the heatmap
ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(limits = c(-2, 2)) # Set the limits for the color scal
######################################################################
#                          FIGURES                                   #
### -------------------------------------------------------------- ###
### HeatMap Plot of Frequency of Plausible 1 - 2 Source Replicates ###
### -------------------------------------------------------------- ###

color_feq_scale_max = max(qpAdm_summary_table_branch_RepSummary_01$mean_plausible_01, 
                          qpAdm_summary_table_aDNA_RepSummary_01$mean_plausible_01,
                          qpAdm_summary_table_branch_RepSummary_05$mean_plausible_05, 
                          qpAdm_summary_table_aDNA_RepSummary_05$mean_plausible_05)
color_feq_scale_min = min(qpAdm_summary_table_branch_RepSummary_01$mean_plausible_01, 
                               qpAdm_summary_table_aDNA_RepSummary_01$mean_plausible_01,
                          qpAdm_summary_table_branch_RepSummary_05$mean_plausible_05, 
                          qpAdm_summary_table_aDNA_RepSummary_05$mean_plausible_05)
  
### -------------------------------------------------------------- ###
### -------------------- 0.01 plausible frequency            ----- ###
### --------- ###
### F2 Branch ###
### --------- ###
Figure_Heatmap__qpAdm_p001_ts_branch <- ggplot(data = qpAdm_summary_table_branch_RepSummary_01,
                                               aes(x=Source1, y=Source2, fill=mean_plausible_01))+
  geom_tile(size=0.5, color = "black")+
  scale_fill_gradient(limits = c(color_feq_scale_min, color_feq_scale_max))+
  geom_text(aes(label=mean_plausible_01), show.legend = F)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1, size=10, face="bold"))+
  theme(axis.text.y = element_text(size=10, face="bold"))+
  ggtitle(paste0("Plausible (p >= 0.01) single and two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2 ", "\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_Heatmap__qpAdm_p001_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_Heatmap__qpAdm_p001_ts_branch, path = out_Fig)


### --------- ###
### sim aDNA  ###
### --------- ###
Figure_Heatmap__qpAdm_p001_aDNA <- ggplot(data = qpAdm_summary_table_aDNA_RepSummary_01,
                                          aes(x=Source1, y=Source2, fill=mean_plausible_01))+
  geom_tile(size=0.5, color = "black")+
  scale_fill_gradient(limits = c(color_feq_scale_min, color_feq_scale_max))+
  geom_text(aes(label=mean_plausible_01), show.legend = F)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1, size=10, face="bold"))+
  theme(axis.text.y = element_text(size=10, face="bold"))+
  ggtitle(paste0("Plausible (p >= 0.01) single and two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA ", "\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_Heatmap__qpAdm_p001_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_Heatmap__qpAdm_p001_aDNA, path = out_Fig)


### -------------------------------------------------------------- ###
### -------------------- 0.05 plausible frequency            ----- ###
### --------- ###
### F2 Branch ###
### --------- ###
Figure_Heatmap__qpAdm_p005_ts_branch <- ggplot(data = qpAdm_summary_table_branch_RepSummary_05,
                                               aes(x=Source1, y=Source2, fill=mean_plausible_05))+
  geom_tile(size=0.5, color = "black")+
  scale_fill_gradient(limits = c(color_feq_scale_min, color_feq_scale_max))+
  geom_text(aes(label=mean_plausible_05), show.legend = F)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1, size=10, face="bold"))+
  theme(axis.text.y = element_text(size=10, face="bold"))+
  ggtitle(paste0("Plausible (p >= 0.05) single and two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2","\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_Heatmap__qpAdm_p005_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_Heatmap__qpAdm_p005_ts_branch, path = out_Fig)

### --------- ###
### sim aDNA  ###
### --------- ###
Figure_Heatmap__qpAdm_p005_aDNA <- ggplot(data = qpAdm_summary_table_aDNA_RepSummary_05,
                                          aes(x=Source1, y=Source2, fill=mean_plausible_05))+
  geom_tile(size=0.5, color = "black")+
  scale_fill_gradient(limits = c(color_feq_scale_min, color_feq_scale_max))+
  geom_text(aes(label=mean_plausible_05), show.legend = F)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1, size=10, face="bold"))+
  theme(axis.text.y = element_text(size=10, face="bold"))+
  ggtitle(paste0("Plausible (p >= 0.05) single and two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA","\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_Heatmap__qpAdm_p005_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_Heatmap__qpAdm_p005_aDNA, path = out_Fig)




#' aDNA and F2 branch together
Figure_Heatmap__qpAdm_p001_p005_aDNA_ts_branch = cowplot::plot_grid(Figure_Heatmap__qpAdm_p001_aDNA, Figure_Heatmap__qpAdm_p001_ts_branch, 
                                                                    Figure_Heatmap__qpAdm_p005_aDNA, Figure_Heatmap__qpAdm_p005_ts_branch, ncol=2)
ggsave(paste0("Figure_Heatmap__qpAdm_p001_p005_aDNA_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_Heatmap__qpAdm_p001_p005_aDNA_ts_branch, path = out_Fig)









####################################################
#                     FIGURES                      #
### ---------------------------------------------###
###          P - VALUES DISTRIBUTION             ###
###                  FIGURES                     ###
### -------------------------------------------- ###
as.factor(branch_qpAdmRotate$Source1)

### ---------------------- ###
###  1 Source Replicates   ###
###                        ###
### -----------------------###
### --------- ###
### F2 Branch ###
### FIGURES   ###
### --------- ###
Figure_qpAdm_1Source_pvalueDist_ts_branch <- branch_qpAdmRotate[numLeftPops == 1] %>% ggplot() +
  geom_point(aes(admix_model_ID, -log10(pvalue),  color=as.character(plausible)))+
  #geom_text_repel(data = branch_qpAdmRotate[numLeftPops == 1], aes(x = admix_model_ID, y = -log10(pvalue), label = branch_qpAdmRotate[numLeftPops == 1, sim]), 
  #                size = 3, vjust = 0, hjust = -0.5)+
  #scale_x_discrete(labels=unique(qpAdm_summary_table_branch[numLeftPops == 1]$Source1))+
  scale_color_manual(name = "Plausible p-value >= 0.01 replicate", values = setNames(c('red', 'black'), c(1,0)))+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), lty=2)+
  theme_bw()+
  coord_cartesian(ylim=c(0,3))+
  xlab("Frequency of plausible qpAdm admixture models")+
  ggtitle(paste0("Plausible (p >= 0.01) one source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2", "\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_qpAdm_1Source_pvalueDist_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_1Source_pvalueDist_ts_branch, path = out_Complx_Fig)

Figure_qpAdm_1Source_pvalueDist_aDNA <- aDNA_qpAdmRotate[numLeftPops == 1] %>% ggplot() +
  geom_point(aes(admix_model_ID, -log10(pvalue),  color=as.character(plausible)))+
  geom_text_repel(data = aDNA_qpAdmRotate[numLeftPops == 1], aes(x = admix_model_ID, y = -log10(pvalue), label = aDNA_qpAdmRotate[numLeftPops == 1]$sim), 
                  size = 3, vjust = 0, hjust = -0.5)+
  #scale_x_discrete(labels=unique(aDNA_qpAdmRotate[numLeftPops == 1]$Source1))+
  scale_color_manual(name = "Plausible p-value >= 0.01 replicate", values = setNames(c('red', 'black'), c(1,0)))+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), lty=2)+
  theme_bw()+
  coord_cartesian(ylim=c(0,3))+
  xlab("Frequency of plausible qpAdm admixture models")+
  ggtitle(paste0("Plausible (p >= 0.01) one source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA", "\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_qpAdm_1Source_pvalueDist_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_1Source_pvalueDist_aDNA, path = out_Complx_Fig)





### ---------------------- ###
###  2 Source Replicates   ###
###                        ###
### -----------------------###

### --------- ###
### F2 Branch ###
### FIGURES   ###
### --------- ###
Figure_qpAdm_2Source_pvalueDist_ts_branch <- ggplot() +
  geom_point(data=branch_qpAdmRotate, aes(admix_model_ID, -log10(pvalue),  color=as.character(plausible)))+
  #geom_text_repel(data = branch_qpAdmRotate, aes(x = admix_model_ID, y = -log10(pvalue), label = branch_qpAdmRotate$sim), 
  #                size = 3, vjust = 0, hjust = -0.5, max.overlaps = Inf)+
  scale_x_discrete(labels=mixedsort(unique(branch_qpAdmRotate$Xlab_01_FullLabel)))+
  scale_color_manual(name = "Plausible p-value >= 0.01 replicate", values = setNames(c('red', 'black'), c('True','False')))+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), lty=2)+
  theme_bw()+
  coord_cartesian(ylim=c(0,3))+
  xlab("Frequency of plausible qpAdm admixture models")+
  ggtitle(paste0("Plausible (p >= 0.01) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2", "\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_qpAdm_2Source_pvalueDist_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_2Source_pvalueDist_ts_branch, path = out_Complx_Fig)
### --------- ###
### sim aDNA  ###
### FIGURES   ###
### --------- ###
Figure_qpAdm_2Source_pvalueDist_aDNA <- ggplot() +
  geom_point(data=qpAdm_summary_table_aDNA_2Source, aes(admix_model_ID, -log10(pvalue),  color=as.character(plausible)))+
  geom_text_repel(data = qpAdm_summary_table_aDNA_2Source, aes(x = admix_model_ID, y = -log10(pvalue), label = qpAdm_summary_table_aDNA_2Source$rep), 
                  size = 3, vjust = 0, hjust = -0.5)+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_aDNA_2Source$Xlab_01_FullLabel)))+
  scale_color_manual(name = "Plausible p-value >= 0.01 replicate", values = setNames(c('red', 'black'), c('True','False')))+
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), lty=2)+
  theme_bw()+
  coord_cartesian(ylim=c(0,3))+
  xlab("Frequency of plausible qpAdm admixture models")+
  ggtitle(paste0("Plausible (p >= 0.01) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA","\n",
                 "Target: ", taret_pop))
ggsave(paste0("Figure_qpAdm_2Source_pvalueDist_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_2Source_pvalueDist_aDNA, path = out_Complx_Fig)


### ------------------ ###
###   FIGURE LEGEND    ###
###                    ###
### ------------------ ###
# Figure__pvalues_aDNA_Legend <- ggplot() +
#   geom_point(data=qpAdm_summary_table_aDNA_2Source, aes(admix_model_ID, -log10(pvalue),  fill=PopLabel))+
#   geom_text_repel(data = qpAdm_summary_table_aDNA_2Source, aes(x = admix_model_ID, y = -log10(pvalue), label = qpAdm_summary_table_aDNA_2Source$rep), 
#                   size = 3, vjust = 0, hjust = -0.5)+
#   theme(legend.key.size = unit(1.5, units = "cm"))+
#   guides(fill=guide_legend(title="qpAdm Model IDs & Source Populations"))
# pvalues_aDNA_Legend <- cowplot::get_legend(Figure__pvalues_aDNA_Legend)
# grid.newpage()
# grid.draw(pvalues_aDNA_Legend)
# ggsave(paste0("pvalues_aDNA_Legend", simModel,"_sim_aDNA",
#               "__","Target_", taret_pop$Pop, "_gen_", taret_pop$Generations,".pdf"), pvalues_aDNA_Legend, path = out_aDNA)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





####################################################
#                  FIGURES                         #
### ---------------------------------------------###
###          ADMXITURE WEIGHTS                   ###
###               FIGURES                        ###
### -------------------------------------------- ###


###            0.05 Significance                   ###
### ---------------------------------------------- ###
###                BRANCH F2                       ###
### ---------------------------------------------- ###
Figure_qpAdm_p005_2Source_weightsDist_ts_branch <- ggplot(data=qpAdm_summary_table_branch_2Source, x=admix_model_ID, y=ifelse(plausible_05 == "True", weight1, NA))+
  geom_boxplot(aes(admix_model_ID, ifelse(plausible_05 == "True", weight1, NA), fill=PopLabel))+
  geom_dotplot(aes(admix_model_ID, ifelse(plausible_05 == "True", weight1, NA)), binaxis='y', stackdir='center', dotsize=0.3)+
  #annotate("point", x = "5", y = 0.9, size=3, colour = "red")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  #geom_text_repel(data = qpAdm_summary_table_branch_2Source[plausible_05 == "True"], aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_branch_2Source[plausible_05 == "True"]$sim),
  #                size = 3, vjust = 0, hjust = -0.5)+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_branch_2Source$Xlab_05_FullLabel)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  ggtitle(paste0("qpAdm weights for plausible (p >= 0.05) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2","\n",
                 "Target: ", taret_pop))+
  theme(legend.key.size = unit(1.5, units = "cm"))
ggsave(paste0("Figure_qpAdm_p005_2Source_weightsDist_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p005_2Source_weightsDist_ts_branch, path = out_Complx_Fig)



### ---------------------------------------------- ###
###                sim aDNA                        ###
### ---------------------------------------------- ###
Figure_qpAdm_p005_2Source_weightsDist_aDNA <- ggplot(data=qpAdm_summary_table_aDNA_2Source, x=admix_model_ID, y=ifelse(plausible_05 == "True", weight1, NA))+
  geom_boxplot(aes(admix_model_ID, ifelse(plausible_05 == "True", weight1, NA), fill=PopLabel))+
  geom_dotplot(aes(admix_model_ID, ifelse(plausible_05 == "True", weight1, NA)), binaxis='y', stackdir='center', dotsize=0.3)+
  #annotate("point", x = "5", y = 0.9, size=3, colour = "red")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  #geom_text_repel(data = qpAdm_summary_table_aDNA_2Source[plausible_05 == "True"], aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_aDNA_2Source[plausible_05 == "True"]$rep),
  #                size = 3, vjust = 0, hjust = -0.5)+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_aDNA_2Source$Xlab_05_FullLabel)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  ggtitle(paste0("qpAdm weights for plausible (p >= 0.05) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA","\n",
                 "Target: ", taret_pop))+
  theme(legend.key.size = unit(1.5, units = "cm"))
ggsave(paste0("Figure_qpAdm_p005_2Source_weightsDist_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p005_2Source_weightsDist_aDNA, path = out_Complx_Fig)


Figure_qpAdm_p005_2Source_weightsDist_aDNA_ts_branch = cowplot::plot_grid(Figure_qpAdm_p005_2Source_weightsDist_aDNA, 
                                                                          Figure_qpAdm_p005_2Source_weightsDist_ts_branch, ncol=1)
ggsave(paste0("Figure_qpAdm_p005_2Source_weightsDist_aDNA_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p005_2Source_weightsDist_aDNA_ts_branch, path = out_Complx_Fig)


###            0.01 Significance                   ###  
### ---------------------------------------------- ###
###                BRANCH F2                       ###
### ---------------------------------------------- ###
Figure_qpAdm_p001_2Source_weightsDist_ts_branch <- ggplot(data=qpAdm_summary_table_branch_2Source, x=admix_model_ID, y=ifelse(plausible == "True", weight1, NA))+
  geom_boxplot(aes(admix_model_ID, ifelse(plausible == "True", weight1, NA), fill=PopLabel))+
  geom_dotplot(aes(admix_model_ID, ifelse(plausible == "True", weight1, NA)), binaxis='y', stackdir='center', dotsize=0.3)+
  #annotate("point", x = "5", y = 0.9, size=3, colour = "red")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  geom_text_repel(data = qpAdm_summary_table_branch_2Source[plausible == "True"], aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_branch_2Source[plausible == "True"]$rep),
                  size = 3, vjust = 0, hjust = -0.5)+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_branch_2Source$Xlab_01_FullLabel)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  ggtitle(paste0("qpAdm weights for plausible (p >= 0.01) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: branch F2","\n",
                 "Target: ", taret_pop))+
  theme(legend.key.size = unit(1.5, units = "cm"))
ggsave(paste0("Figure_qpAdm_p001_2Source_weightsDist_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p001_2Source_weightsDist_ts_branch, path = out_Complx_Fig)



### ---------------------------------------------- ###
###                sim aDNA                        ###
### ---------------------------------------------- ###
Figure_qpAdm_p001_2Source_weightsDist_aDNA <- ggplot(data=qpAdm_summary_table_aDNA_2Source, x=admix_model_ID, y=ifelse(plausible == "True", weight1, NA))+
  geom_boxplot(aes(admix_model_ID, ifelse(plausible == "True", weight1, NA), fill=PopLabel))+
  geom_dotplot(aes(admix_model_ID, ifelse(plausible == "True", weight1, NA)), binaxis='y', stackdir='center', dotsize=0.3)+
  #annotate("point", x = "5", y = 0.9, size=3, colour = "red")+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  geom_text_repel(data = qpAdm_summary_table_aDNA_2Source[plausible == "True"], aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_aDNA_2Source[plausible == "True"]$rep),
                  size = 3, vjust = 0, hjust = -0.5)+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_aDNA_2Source$Xlab_01_FullLabel)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  ggtitle(paste0("qpAdm weights for plausible (p >= 0.01) two source qpAdm replicates ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Data type: sim aDNA","\n",
                 "Target: ", taret_pop))+
  theme(legend.key.size = unit(1.5, units = "cm"))
ggsave(paste0("Figure_qpAdm_p001_2Source_weightsDist_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p001_2Source_weightsDist_ts_branch, path = out_Complx_Fig)


#' aDNA and F2 branch together
Figure_qpAdm_p001_2Source_weightsDist_aDNA_ts_branch = cowplot::plot_grid(Figure_qpAdm_p001_2Source_weightsDist_aDNA, 
                                                                          Figure_qpAdm_p001_2Source_weightsDist_ts_branch, ncol=1)
ggsave(paste0("Figure_qpAdm_p001_2Source_weightsDist_aDNA_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 15, height = 10, units = "in", dpi=300, Figure_qpAdm_p001_2Source_weightsDist_aDNA_ts_branch, path = out_Complx_Fig)


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





####################################################
#                  FIGURES                         #
### ---------------------------------------------###
###          ADMXITURE WEIGHTS & p-values        ###
###               FIGURES                        ###
### -------------------------------------------- ###
Figure_qpAdm_2Source_weights_pvalue_Dist_aDNA <- ggplot(data=qpAdm_summary_table_aDNA_2Source, aes(x=admix_model_ID, y=weight1))+
  geom_violin()+
  coord_cartesian(ylim=c(-10,10))+
  geom_hline(yintercept = c(1, 0), linetype = "dashed")+
  geom_jitter(data=qpAdm_summary_table_aDNA_2Source, aes(shape = pvalue >= 0.05, color = plausible), size = 3, width=0.15, alpha=0.5)+
  #geom_text_repel(data = qpAdm_summary_table_aDNA_2Source, aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_aDNA_2Source$rep),
  #                size = 3, vjust = 0, hjust = -0.5)+
  theme_bw()+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_aDNA_2Source$Xlab_fq_pWeights)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  theme(plot.title = element_text(size = 10, face = "bold"))+
  labs(subtitle = paste0("Data type: sim aDNA", "\n",
                         "Mean all-sources QTP = ", mean(qpAdm_summary_table_aDNA[admix_model_ID == 1]$QTP_p005_weight01)))
ggsave(paste0("Figure_qpAdm_2Source_weights_pvalue_Dist_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 20, height = 10, units = "in", dpi=300, Figure_qpAdm_2Source_weights_pvalue_Dist_aDNA, path = out_Complx_Fig)


Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch <- ggplot(data=qpAdm_summary_table_branch_2Source, aes(x=admix_model_ID, y=weight1))+
  geom_violin()+
  coord_cartesian(ylim=c(-10,10))+
  geom_hline(yintercept = c(1, 0), linetype = "dashed")+
  geom_jitter(data=qpAdm_summary_table_branch_2Source, aes(shape = pvalue >= 0.05, color = plausible), size = 3, width=0.15, alpha=0.5)+
  geom_text_repel(data = qpAdm_summary_table_branch_2Source, aes(x = admix_model_ID, y = weight1, label = qpAdm_summary_table_branch_2Source$rep),
                  size = 3, vjust = 0, hjust = -0.5)+
  theme_bw()+
  scale_x_discrete(labels=mixedsort(unique(qpAdm_summary_table_branch_2Source$Xlab_fq_pWeights)))+
  ylab("P1 weight estimate")+
  xlab("qpAdm admixture models")+
  ggtitle(paste0("qpAdm replicate two source model weights ", "\n" , 
                 "Simulation model: ", simModel, "\n", 
                 "Target: ", taret_pop,"\n"))+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))+
  labs(subtitle = paste0("Data type: branch F2", "\n",
                         "Mean all-sources QTP = ", mean(qpAdm_summary_table_branch[admix_model_ID == 1]$QTP_p005_weight01)))
ggsave(paste0("Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 20, height = 10, units = "in", dpi=300, Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch, path = out_Complx_Fig)


Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch_aDNA = cowplot::plot_grid(Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch, 
                                                                             Figure_qpAdm_2Source_weights_pvalue_Dist_aDNA, ncol=1)
ggsave(paste0("Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch_aDNA", "__", simModel, "__" , taret_pop, ".pdf"),  
       width = 20, height = 10, units = "in", dpi=300, Figure_qpAdm_2Source_weights_pvalue_Dist_ts_branch_aDNA, path = out_Complx_Fig)


