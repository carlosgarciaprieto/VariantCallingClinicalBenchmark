#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure S2 
# Usage: change data directory (dataDir) variable and execute the script step-by-step
################################################################################
### Libraries
library(maftools)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(caret)
library(tidyr)
library(tidytext)
library(viridis)
library(ggthemes)
library(season)
library(gridExtra)
library(scales)
library(ggalluvial)
library(ComplexUpset)
library(dplyr)
library(readxl)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(devtools)

################################################################################
## Set working and data directory
dataDir<-"/Users/carlosgarciaprieto/Bioinformatics"
setwd(dataDir)
#===============================================================================
### MAIN.
################################################################################
### Step 1: Read intOGen two-callers intersection results
################################################################################

## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          "LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          "THCA","THYM","UCEC","UCS","UVM")

caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")

## Read intOGen two-caller intersection results saved under Intogen_two_caller folder
# We selected as cancer driver genes those with a significant q.value < 0.5
path_two<-paste(dataDir,"Intogen_two_callers",sep = "/")
files_two<-list.files(path=path_two, pattern = "*05.out", full.names = TRUE)
file_list_two <- sapply(files_two, fread, simplify=FALSE)
# We defined as cancer driver genes those within TIER1 & TIER2 
genes_list_two<-lapply(file_list_two, function(x) {subset(x[[1]],x[[2]]<3)})
genes_two<-reshape2::melt(genes_list_two)
genes_two$L1<-gsub(paste0(print(path_two),"/"),"",genes_two$L1)
genes_two[,2]<-sapply(strsplit(genes_two[,2], "\\."), `[`, 1)
colnames(genes_two)<-c("SYMBOL","Cohort")
genes_two$Cancer<-sapply(strsplit(genes_two[,2], "\\_"), `[`, 2)
genes_two$Caller<-sapply(strsplit(genes_two[,2], "\\_"), `[`, 3)
# We set a key value to reshape data frame
genes_two$Key<-1
genes_two_dcast<-reshape2::dcast(genes_two, SYMBOL + Cancer ~ Caller, value.var = "Key")
genes_two_dcast[is.na(genes_two_dcast)] <- 0
genes_two_dcast$Sum<-apply(genes_two_dcast[,3:9],1,function(x) sum(x))
head(genes_two_dcast)
summary(genes_two_dcast$Sum)
# Check number of cancer driver genes per cancer and caller
genes_two_aggregate<-aggregate(Key ~ Caller + Cancer, data = genes_two, FUN = sum)

################################################################################
### Step 2: Performance metrics of two-callers strategies to detect cancer driver genes  
################################################################################

# Format and save data
genes_two_dcast$SYMBOL_CANCER_TYPE<-paste(genes_two_dcast$SYMBOL,genes_two_dcast$Cancer,sep = "_")
saveRDS(genes_two_dcast,paste(dataDir,"Resources/genes_two_dcast.rds",sep = "/"))
# Read data
genes_two_dcast<-readRDS(paste(dataDir,"Resources/genes_two_dcast.rds",sep = "/"))
intogen<-readRDS(paste(dataDir,"Resources/intogen_TCGA.rds",sep = "/"))
mc3<-readRDS(paste(dataDir,"Resources/mc3.rds",sep = "/"))
head(genes_two_dcast)
head(intogen)
head(mc3)
# Select cancer driver genes in ACC,BLCA,BRCA,PRAD and UCEC
cancers_two<-c("ACC","BLCA","BRCA","PRAD","UCEC")
intogen_two<-intogen[intogen$CANCER_TYPE %in% cancers_two,]
mc3_two<-mc3[mc3$CANCER_TYPE %in% cancers_two,]
# Merge data at gene and cancer level
metrics_two<-Reduce(function(x,y) merge(x,y,by="SYMBOL_CANCER_TYPE",all=TRUE) ,list(genes_two_dcast[,c(3:9,11)],intogen_two[,-c(1:2)],mc3_two[,-c(1:2)]))
metrics_two[is.na(metrics_two)] <- 0
metrics_two<-metrics_two[,c(1:8,12,15,9:11,13:14)]
metrics_two$Total<-apply(metrics_two[,2:8],1,function(x) sum(x))
metrics_two<-metrics_two[,c(1:8,16,9:15)]
metrics_two$SYMBOL<-sapply(strsplit(metrics_two[,1], "\\_"), `[`, 1)
metrics_two$CANCER<-sapply(strsplit(metrics_two[,1], "\\_"), `[`, 2)
metrics_two<-metrics_two[,c(17:18,1:16)]
head(metrics_two)
callers_two<-c("Consensus2","MuSE&MuTect2","MuSE&SomaticSniper","MuSE&VarScan2","MuTect2&SomaticSniper","MuTect2&VarScan2","SomaticSniper&VarScan2")
colnames(metrics_two)[4:10]<-callers_two
# Save results
saveRDS(metrics_two,paste(dataDir,"Resources/metrics_two_TCGA.rds",sep = "/"))
#Write report
writexl::write_xlsx(metrics_two,paste(dataDir,"Supplementary_files/Supplementary_file_3_intOGen_results_combinations_of_two.xlsx",sep = "/"))

## Compute metrics per cancer type
# Read data
metrics_two<-readRDS(paste(dataDir,"Resources/metrics_two_TCGA.rds",sep = "/"))
head(metrics_two)
# Prepare data
callers_two<-c("Consensus2","MuSE&MuTect2","MuSE&SomaticSniper","MuSE&VarScan2","MuTect2&SomaticSniper","MuTect2&VarScan2","SomaticSniper&VarScan2")
# Metrics per cancer type
confusion_matrix_metrics_two_append<-NULL
for(i in 1:length(cancers_two)){
  print(i)
  print(cancers_two[i])
  for(k in 1:length(callers_two)){
    print(paste(i,cancers_two[i],callers_two[k],sep = "_"))
    variant_callers_two<-NULL
    cancers_two_type<-NULL
    confusion_matrix_intogen<-NULL
    confusion_matrix_mc3<-NULL
    confusion_matrix_metrics_two<-NULL
    variant_callers_two<-callers_two[k]
    cancers_two_type<-cancers_two[i]
    confusion_matrix_intogen<-confusionMatrix(data=factor(metrics_two[,variant_callers_two][metrics_two$CANCER==cancers_two_type]),reference=factor(metrics_two$intOGen[metrics_two$CANCER==cancers_two_type]), mode = "prec_recall", positive = "1")
    confusion_matrix_mc3<-confusionMatrix(data=factor(metrics_two[,variant_callers_two][metrics_two$CANCER==cancers_two_type]),reference=factor(metrics_two$mc3[metrics_two$CANCER==cancers_two_type]), mode = "prec_recall", positive = "1")
    confusion_matrix_metrics_two<-as.data.frame(rbind(confusion_matrix_intogen$byClass[5:7],confusion_matrix_mc3$byClass[5:7]))
    confusion_matrix_metrics_two$Caller<-callers_two[k]
    confusion_matrix_metrics_two$Cancer<-cancers_two[i]
    confusion_matrix_metrics_two$Truthset<-c("intOGen","MC3")
    confusion_matrix_metrics_two_append<-rbind(confusion_matrix_metrics_two_append,confusion_matrix_metrics_two)
  }
}

################################################################################
### Step 3: Plot Figure S2
################################################################################

## Plot metrics per cancer type
# Select reference set
confusion_matrix_metrics_two_append_intogen<-confusion_matrix_metrics_two_append[confusion_matrix_metrics_two_append$Truthset=="intOGen",]
confusion_matrix_metrics_two_append_mc3<-confusion_matrix_metrics_two_append[confusion_matrix_metrics_two_append$Truthset=="MC3",]
head(confusion_matrix_metrics_two_append_intogen)
head(confusion_matrix_metrics_two_append_mc3)
summary(confusion_matrix_metrics_two_append_intogen)
summary(confusion_matrix_metrics_two_append_mc3)

## Plot metrics per cancer type as boxplots
# intOGen reference set
metrics_two_boxplot_intogen <- reshape2::melt(confusion_matrix_metrics_two_append_intogen, id.vars = c("Caller","Truthset","Cancer"),variable.name = "Metrics", value.name = "Score")
metrics_two_boxplot_intogen$Metrics<-gsub("F1","F1-Score",metrics_two_boxplot_intogen$Metrics)
metrics_two_boxplot_intogen$Metrics<-factor(metrics_two_boxplot_intogen$Metrics, levels = c("Precision","Recall","F1-Score"))
head(metrics_two_boxplot_intogen)
# Plot Figure S2 (top)
m_two_boxplot_intogen<-ggplot(data = metrics_two_boxplot_intogen, aes(reorder_within(Caller,Score,Metrics,mean),Score,)) + theme_bw(base_size = 12) + geom_boxplot(aes(fill = Caller),width = 0.8, outlier.colour=NA) + scale_x_reordered() + ylab ("Score") + xlab (element_blank()) + facet_grid (as.factor(metrics_two_boxplot_intogen$Truthset)~as.factor (metrics_two_boxplot_intogen$Metrics), scales = "free_x") + theme(axis.text.x = element_text (angle = 45, hjust = 1)) + theme(aspect.ratio = 1, legend.position = "bottom", legend.margin=margin(t=-15),axis.ticks.x=element_blank(), axis.text.x = element_blank(),legend.title = element_blank(),legend.key.size = unit(1,"line"), legend.text = element_text(size = 10, colour = "black", face = "bold")) + guides(fill = guide_legend(nrow = 3)) + scale_fill_manual(values=c("#8F7700FF","chocolate","darkblue","darkcyan","lightgoldenrod","azure3","cyan")) + scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + theme(strip.background = element_rect(colour="black", fill="antiquewhite",size=1, linetype="solid")) + theme(axis.text.y = element_text(color="black", size=12, face="bold"),axis.title.y = element_text(color="black", size=14, face="bold"),strip.text.x = element_text(size=12, color="black",face="bold"),strip.text.y = element_text(size=12, color="black",face="bold")) 
ggsave(plot = m_two_boxplot_intogen, width = 7, height = 7, dpi = 300,filename = paste(dataDir,"Figures/FigureS2_intogen.png",sep = "/"))
saveRDS(m_two_boxplot_intogen,paste(dataDir,"Figures/FigureS2_intogen.rds",sep = "/"))
# mc3 reference set
metrics_two_boxplot_mc3 <- reshape2::melt(confusion_matrix_metrics_two_append_mc3, id.vars = c("Caller","Truthset","Cancer"),variable.name = "Metrics", value.name = "Score")
metrics_two_boxplot_mc3$Metrics<-gsub("F1","F1-Score",metrics_two_boxplot_mc3$Metrics)
metrics_two_boxplot_mc3$Metrics<-factor(metrics_two_boxplot_mc3$Metrics, levels = c("Precision","Recall","F1-Score"))
head(metrics_two_boxplot_mc3)
# Plot Figure S2 (bottom)
m_two_boxplot_mc3<-ggplot(data = metrics_two_boxplot_mc3, aes(reorder_within(Caller,Score,Metrics,mean),Score,)) + theme_bw(base_size = 12) + geom_boxplot(aes(fill = Caller),width = 0.8, outlier.colour=NA) + scale_x_reordered() + ylab ("Score") + xlab (element_blank()) + facet_grid (as.factor(metrics_two_boxplot_mc3$Truthset)~as.factor (metrics_two_boxplot_mc3$Metrics), scales = "free_x") + theme(axis.text.x = element_text (angle = 45, hjust = 1)) + theme(aspect.ratio = 1, legend.position = "bottom", legend.margin=margin(t=-15),axis.ticks.x=element_blank(), axis.text.x = element_blank(),legend.title = element_blank(),legend.key.size = unit(1,"line"), legend.text = element_text(size = 10, colour = "black", face = "bold")) + guides(fill = guide_legend(nrow = 3)) + scale_fill_manual(values=c("#8F7700FF","chocolate","darkblue","darkcyan","lightgoldenrod","azure3","cyan")) + scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + theme(strip.background = element_rect(colour="black", fill="antiquewhite",size=1, linetype="solid")) + theme(axis.text.y = element_text(color="black", size=12, face="bold"),axis.title.y = element_text(color="black", size=14, face="bold"),strip.text.x = element_text(size=12, color="black",face="bold"),strip.text.y = element_text(size=12, color="black",face="bold")) 
ggsave(plot = m_two_boxplot_mc3, width = 7, height = 7, dpi = 300,filename = paste(dataDir,"Figures/FigureS2_mc3.png",sep = "/"))
saveRDS(m_two_boxplot_mc3,paste(dataDir,"Figures/FigureS2_mc3.rds",sep = "/"))
#===============================================================================