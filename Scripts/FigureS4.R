#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Supplementary Figure S4
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
### Step 1: Read performance metrics to detect cancer driver genes with intOGen
################################################################################

# Read metrics results
metrics<-readRDS(paste(dataDir,"Resources/metrics_TCGA.rds",sep = "/"))
head(metrics)
# Format group 
metrics$Group<-ifelse(metrics$CANCER %in% Hematologic_lymphatic,"Hematologic and lymphatic",ifelse(metrics$CANCER %in% Gynecologic, "Gynecologic", ifelse(metrics$CANCER %in% Urologic, "Urologic", ifelse(metrics$CANCER %in% Endocrine, "Endocrine",ifelse(metrics$CANCER %in% Gastrointestinal,"Gastrointestinal",ifelse(metrics$CANCER %in% Thoracic,"Thoracic",ifelse(metrics$CANCER %in% CNS,"Central Nervous System","Other")))))))
#Transform blank spaces to NA!!!
metrics[metrics==""] <- NA
head(metrics)
#Format data
colnames(metrics)[c(2,4:10,13)]<-c("Cancer","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2","MC3")
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers_driver_genes<-c("MC3","intOGen",callers)
metrics[c(callers_driver_genes)]<-metrics[c(callers_driver_genes)] != 0
head(metrics)
groups<-sort(unique(metrics$Group))

################################################################################
### Step 2: Plot Figure S4
################################################################################

# Select color blind friendly palette
cbPalette <- c("#E69F00","#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Plot Figure S4 (part 1)
Figure_intogen_results1<-
  (
    upset(metrics[metrics$Group==groups[1],], callers_driver_genes, name = NULL ,width_ratio=0.2, height_ratio = 0.8, n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('GBM'=cbPalette[1], 'LGG'=cbPalette[2]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=130)+theme(legend.position = "none")+scale_fill_manual(values=c('GBM'=cbPalette[1], 'LGG'=cbPalette[2]))))
    + ggtitle("Central Nervous System") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[2],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8, n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('ACC'=cbPalette[1], 'THCA'=cbPalette[2]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=130)+theme(legend.position = "none")+scale_fill_manual(values=c('ACC'=cbPalette[1], 'THCA'=cbPalette[2]))))
    + ggtitle(as.character(groups[2])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[3],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8,  n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('CHOL'=cbPalette[1], 'COAD'=cbPalette[2],'ESCA'=cbPalette[3],'LIHC'=cbPalette[4], 'PAAD'=cbPalette[5], 'READ'=cbPalette[6], 'STAD'=cbPalette[7]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.1)+ expand_limits(y=380)+theme(legend.position = "none")+scale_fill_manual(values=c('CHOL'=cbPalette[1], 'COAD'=cbPalette[2],'ESCA'=cbPalette[3],'LIHC'=cbPalette[4], 'PAAD'=cbPalette[5], 'READ'=cbPalette[6], 'STAD'=cbPalette[7]))))
    + ggtitle(as.character(groups[3])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[4],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8,  n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('BRCA'=cbPalette[1], 'CESC'=cbPalette[2],'OV'=cbPalette[3],'UCEC'=cbPalette[4]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=265)+theme(legend.position = "none")+scale_fill_manual(values=c('BRCA'=cbPalette[1], 'CESC'=cbPalette[2],'OV'=cbPalette[3],'UCEC'=cbPalette[4]))))
    + ggtitle(as.character(groups[4])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_intogen_results1, width = 11, height = 8, dpi = 300, filename = paste(dataDir,"Figures/FigureS4_1.png",sep = "/"))
#Save data
saveRDS(Figure_intogen_results1,paste(dataDir,"Figures/FigureS4_1.rds",sep = "/"))
# Plot Figure S4 (part 2)
Figure_intogen_results2<-
  (
    upset(metrics[metrics$Group==groups[5],], callers_driver_genes, name = NULL ,width_ratio=0.2, height_ratio = 0.8, n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('DLBC'=cbPalette[1], 'LAML'=cbPalette[2], 'THYM'=cbPalette[3]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=175)+theme(legend.position = "none")+scale_fill_manual(values=c('DLBC'=cbPalette[1], 'LAML'=cbPalette[2], 'THYM'=cbPalette[3]))))
    + ggtitle(as.character(groups[5])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[7],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8,  n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('LUAD'=cbPalette[1], 'LUSC'=cbPalette[2],'MESO'=cbPalette[3]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=150)+theme(legend.position = "none")+scale_fill_manual(values=c('LUAD'=cbPalette[1], 'LUSC'=cbPalette[2],'MESO'=cbPalette[3]))))
    + ggtitle(as.character(groups[7])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[8],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8, n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('BLCA'=cbPalette[1], 'KICH'=cbPalette[2], 'KIRC'=cbPalette[3], 'KIRP'=cbPalette[4], 'PRAD'=cbPalette[5], 'TGCT'=cbPalette[6]))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=325)+theme(legend.position = "none")+scale_fill_manual(values=c('BLCA'=cbPalette[1], 'KICH'=cbPalette[2], 'KIRC'=cbPalette[3], 'KIRP'=cbPalette[4], 'PRAD'=cbPalette[5], 'TGCT'=cbPalette[6]))))
    + ggtitle(as.character(groups[8])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(metrics[metrics$Group==groups[6],], callers_driver_genes, name=NULL, width_ratio=0.2, height_ratio = 0.8,  n_intersections = 15, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('HNSC'=cbPalette[1], 'PCPG'=cbPalette[2],'SARC'=cbPalette[3],'SKCM'=cbPalette[4], 'UCS'=cbPalette[5], 'UVM'=cbPalette[6])) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.25)+ expand_limits(y=245)+theme(legend.position = "none")+scale_fill_manual(values=c('HNSC'=cbPalette[1], 'PCPG'=cbPalette[2],'SARC'=cbPalette[3],'SKCM'=cbPalette[4], 'UCS'=cbPalette[5], 'UVM'=cbPalette[6]))))
    + ggtitle(as.character(groups[6])) + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_intogen_results2, width = 11, height = 8, dpi = 300, filename = paste(dataDir,"Figures/FigureS4_2.png",sep = "/"))
# Save data
saveRDS(Figure_intogen_results2,paste(dataDir,"Figures/FigureS4_2.rds",sep = "/"))
#===============================================================================