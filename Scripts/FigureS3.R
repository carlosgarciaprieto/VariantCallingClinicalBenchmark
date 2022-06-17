#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure S3
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
### Step 1: Plot best performing variant calling strategy (Figure S3) 
################################################################################

## Classify cancer types
Hematologic_lymphatic<-c("LAML","DLBC","THYM")
Gynecologic<-c("OV","UCEC","CESC","BRCA")
Urologic<-c("BLCA","PRAD","TGCT","KIRC","KICH","KIRP")
Endocrine<-c("THCA","ACC")
Core_gastrointestinal<-c("ESCA","STAD","COAD","READ")
Developmental_gastrointestinal<-c("LIHC","PAAD","CHOL")
Gastrointestinal<-c(Core_gastrointestinal,Developmental_gastrointestinal)
Head_and_neck<-c("HNSC")
Thoracic<-c("LUAD","LUSC","MESO")
CNS<-c("GBM","LGG")
Soft_tissue<-c("SARC","UCS")
Neural_crest_derived_tissues<-c("PCPG")
Skin<-"SKCM"
Eye<-"UVM"

## Alluvial plot (Figure S3)
# Read data
alluvial<-readRDS(paste(dataDir,"Resources/alluvial.rds",sep = "/"))
# Prepare data
alluvial$Group<-ifelse(alluvial$Cancer %in% Hematologic_lymphatic,"Hematologic and lymphatic",ifelse(alluvial$Cancer %in% Gynecologic, "Gynecologic", ifelse(alluvial$Cancer %in% Urologic, "Urologic", ifelse(alluvial$Cancer %in% Endocrine, "Endocrine",ifelse(alluvial$Cancer %in% Gastrointestinal,"Gastrointestinal",ifelse(alluvial$Cancer %in% Thoracic,"Thoracic",ifelse(alluvial$Cancer %in% CNS,"Central Nervous System",ifelse(alluvial$Cancer %in% Soft_tissue,"Soft tissue","Other"))))))))
alluvial$Group<-factor(alluvial$Group,levels = c("Central Nervous System","Endocrine","Gastrointestinal","Gynecologic","Hematologic and lymphatic","Soft tissue","Thoracic","Urologic","Other"))
alluvial_intogen<-alluvial[alluvial$Truthset=="intOGen",]
alluvial_mc3<-alluvial[alluvial$Truthset=="MC3",]

## mc3 reference set 
# Precision
alluvial_mc3_plot_precision<-ggplot(data = alluvial_mc3,
                      aes(axis2 = Group, axis1 = `Precision`)) +
                      scale_x_discrete(limits = c("Precision","Group"), expand = c(.0, .0)) +
                      geom_alluvium(aes(fill = `Precision`)) +
                      geom_stratum(color = "black") +
                      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                      theme_classic() +
                      scale_fill_manual(values = c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                      ggtitle("Best Caller")
ggsave(plot = alluvial_mc3_plot_precision, width = 15.4, height = 5.5, dpi = 300, filename = paste(dataDir,"Figures/FigureS3_precision_mc3.png",sep = "/"))
# Recall
alluvial_mc3_plot_recall<-ggplot(data = alluvial_mc3,
                      aes(axis2 = Group, axis1 = `Recall`)) +
                      scale_x_discrete(limits = c("Recall","Group"), expand = c(.0, .0)) +
                      geom_alluvium(aes(fill = `Recall`)) +
                      geom_stratum(color = "black") +
                      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                      theme_classic() +
                      scale_fill_manual(values = c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                      ggtitle("Best Caller")
ggsave(plot = alluvial_mc3_plot_recall, width = 15.4, height = 5.5, dpi = 300, filename = paste(dataDir,"Figures/FigureS3_recall_mc3.png",sep = "/"))

## intOGen reference set
# Precision
alluvial_intogen_plot_precision<-ggplot(data = alluvial_intogen,
                              aes(axis2 = Group, axis1 = `Precision`)) +
                              scale_x_discrete(limits = c("Precision","Group"), expand = c(.0, .0)) +
                              geom_alluvium(aes(fill = `Precision`)) +
                              geom_stratum(color = "black") +
                              geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                              theme_classic() +
                              scale_fill_manual(values = c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                              ggtitle("Best Caller intOGen")
ggsave(plot = alluvial_intogen_plot_precision, width = 9.5, height = 5, dpi = 300, filename = paste(dataDir,"Figures/FigureS3_precision_intogen.png",sep = "/"))
# Recall
alluvial_intogen_plot_recall<-ggplot(data = alluvial_intogen,
                                        aes(axis2 = Group, axis1 = `Recall`)) +
                                        scale_x_discrete(limits = c("Recall","Group"), expand = c(.0, .0)) +
                                        geom_alluvium(aes(fill = `Recall`)) +
                                        geom_stratum(color = "black") +
                                        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                                        theme_classic() +
                                        scale_fill_manual(values = c("#8F7700FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                                        ggtitle("Best Caller intOGen")
ggsave(plot = alluvial_intogen_plot_recall, width = 9.5, height = 5, dpi = 300, filename = paste(dataDir,"Figures/FigureS3_precision_intogen.png",sep = "/"))

# Save results
saveRDS(alluvial_mc3_plot_precision,paste(dataDir,"Figures_Paper/FigureS3_precision_mc3.rds",sep = "/"))
saveRDS(alluvial_mc3_plot_recall,paste(dataDir,"Figures_Paper/FigureS3_recall_mc3.rds",sep = "/"))
saveRDS(alluvial_intogen_plot_precision,paste(dataDir,"Figures_Paper/FigureS3_precision_intogen.rds",sep = "/"))
saveRDS(alluvial_intogen_plot_recall,paste(dataDir,"Figures_Paper/FigureS3_recall_intogen.rds",sep = "/"))
#===============================================================================