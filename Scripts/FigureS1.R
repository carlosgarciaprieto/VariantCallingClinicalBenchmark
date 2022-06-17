#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Supplementary Figure S1
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
### Step 1: Read intOGen downsampling experiment results
################################################################################

# We define as cancer driver genes those with a significant q.value < 0.5
path<-paste(dataDir,"Intogen_Downsampling",sep = "/")
files<-list.files(path=path, pattern = "*05.out", full.names = TRUE)
file_list <- sapply(files, fread, simplify=FALSE)
# We select TIER1 & TIER2 cancer driver genes
downsampling_genes_list<-lapply(file_list, function(x) {subset(x[[1]],x[[2]]<3)})
# We reshape data
downsampling_genes<-reshape2::melt(downsampling_genes_list)
downsampling_genes$L1<-gsub(paste0(print(path),"/"),"",downsampling_genes$L1)
downsampling_genes[,2]<-sapply(strsplit(downsampling_genes[,2], "\\."), `[`, 1)
colnames(downsampling_genes)<-c("SYMBOL","Cohort")
downsampling_genes$Cancer<-sapply(strsplit(downsampling_genes[,2], "\\_"), `[`, 1)
downsampling_genes$Caller<-sapply(strsplit(downsampling_genes[,2], "\\_"), `[`, 2)
downsampling_genes$Percentage<-sapply(strsplit(downsampling_genes[,2], "\\_"), `[`, 3)
downsampling_genes$Iteration<-sapply(strsplit(downsampling_genes[,2], "\\_"), `[`, 4)
# We set a key value to reshape data frame
downsampling_genes$Key<-1

## We add the number of BRCA samples corresponding to each downsampling experiment (25%, 50% and 75%)
# 251 samples (25%)
# 496 samples (50%)
# 745 samples (75)
downsampling_genes$Samples<-ifelse(downsampling_genes$Percentage==25,251,ifelse(downsampling_genes$Percentage==50,496,745))
table(downsampling_genes$Percentage,downsampling_genes$Samples)
# Reshape data
downsampling_genes_dcast<-reshape2::dcast(downsampling_genes, SYMBOL + Cancer  + Percentage + Iteration + Samples ~ Caller, value.var = "Key")
downsampling_genes_dcast[is.na(downsampling_genes_dcast)] <- 0
downsampling_genes_dcast$Sum<-apply(downsampling_genes_dcast[,6:12],1,function(x) sum(x))
head(downsampling_genes_dcast)
summary(downsampling_genes_dcast$Sum)
# Check number of cancer driver genes per cancer and caller
downsampling_genes_aggregate<-aggregate(Key ~ Caller + Cancer + Percentage + Iteration + Samples, data = downsampling_genes, FUN = sum)

################################################################################
### Step 2: Prepare data for plot
################################################################################

# Read data
genes_burden<-readRDS(paste(dataDir,"Resources/genes_burden_size.rds",sep = "/"))
head(genes_burden)
head(downsampling_genes_aggregate)
# Prepare data
genes_burden_BRCA<-genes_burden[genes_burden$Cancer == "BRCA",]
genes_burden_BRCA$Percentage<-100
genes_burden_BRCA$Iteration<-1
genes_burden_BRCA<-genes_burden_BRCA[,c(3,1,9,10,7,4)]
head(genes_burden_BRCA)
downsampling_genes_aggregate_total<-rbind(downsampling_genes_aggregate,genes_burden_BRCA)
tail(downsampling_genes_aggregate_total)
# Save data
saveRDS(downsampling_genes_aggregate_total,paste(dataDir,"downsampling_genes_aggregate_total.rds",sep = "/"))

################################################################################
### Step 3: Plot Figure S1
################################################################################

# Read data
downsampling_genes_aggregate_total<-readRDS(paste(dataDir,"Resources/downsampling_genes_aggregate_total.rds",sep = "/"))
max(downsampling_genes_aggregate_total$Key)
# Prepare data
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "CONSENSUS2"]<-"Consensus2"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "CONSENSUS3"]<-"Consensus3"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "MUSE"]<-"MuSE"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "MUTECT2"]<-"MuTect2"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "SOMATICSNIPER"]<-"SomaticSniper"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "UNION"]<-"Union"
downsampling_genes_aggregate_total$Caller[downsampling_genes_aggregate_total$Caller == "VARSCAN2"]<-"VarScan2"
# Plot Figure S1
set.seed(24)
FigureDownsampling<-ggplot(data = downsampling_genes_aggregate_total, aes(y = Key, x = Samples)) + theme_bw(base_size = 12) + geom_point(aes(colour = Caller, shape = Iteration), size = 3, alpha = 0.85, position=position_jitter(h=0, w=20))  + ylab ("intOGen detected driver genes") + xlab ("Number of samples") + stat_cor(method = "spearman", alternative = "two.sided", cor.coef.name = "rho", size = 5) + stat_smooth(method = "lm", colour = "black", size=0.5, se = T) + scale_color_manual(values=c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) + theme(aspect.ratio = 1, legend.position = "right") + scale_shape_manual(values=c(19,15,17)) + theme(axis.title = element_text(size = 14, face = "bold"),axis.text = element_text(size = 12), legend.text = element_text(size = 11))
ggsave(plot = FigureDownsampling, width = 6.5, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/FigureS1.png",sep = "/"))
# Save data
saveRDS(FigureDownsampling,paste(dataDir,"Figures/FigureS1.rds",sep = "/"))
#===============================================================================