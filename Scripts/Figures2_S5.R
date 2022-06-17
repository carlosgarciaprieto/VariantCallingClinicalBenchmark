#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure 2 and Figure S5
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
### Step 1: Read intOGen results
################################################################################

## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          "LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          "THCA","THYM","UCEC","UCS","UVM")

caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")

## Read intOGen results saved under Intogen folder
# We selected as cancer driver genes those with a significant q.value < 0.5
path<-paste(dataDir,"Intogen",sep = "/")
files<-list.files(path=path, pattern = "*05.out", full.names = TRUE)
file_list <- sapply(files, fread, simplify=FALSE)
# We defined as cancer driver genes those within TIER1 & TIER2 
genes_list<-lapply(file_list, function(x) {subset(x[[1]],x[[2]]<3)})
genes<-reshape2::melt(genes_list)
genes$L1<-gsub(paste0(print(path),"/"),"",genes$L1)
genes[,2]<-sapply(strsplit(genes[,2], "\\."), `[`, 1)
colnames(genes)<-c("SYMBOL","Cohort")
genes$Cancer<-sapply(strsplit(genes[,2], "\\_"), `[`, 2)
genes$Caller<-sapply(strsplit(genes[,2], "\\_"), `[`, 3)
# We set a key value to reshape data frame
genes$Key<-1
genes_dcast<-reshape2::dcast(genes, SYMBOL + Cancer ~ Caller, value.var = "Key")
genes_dcast[is.na(genes_dcast)] <- 0
genes_dcast$Sum<-apply(genes_dcast[,3:9],1,function(x) sum(x))
head(genes_dcast)
summary(genes_dcast$Sum)
# Check number of driver genes per cancer and caller
genes_aggregate<-aggregate(Key ~ Caller + Cancer, data = genes, FUN = sum)
View(genes_aggregate)

################################################################################
### Step 2: Compute Tumor Mutational Burden (TMB)
################################################################################

## COMPUTE TMB 
# For each of the 33 TCGA cancer types, we downloaded the four different Somatic aggregated MAF files with all the somatic mutations for each variant caller (MuSE, MuTect2, SomaticSniper and VarScan2), created the Consensus of 2, Consensus of 3 and Union files, and saved them under Maf_files folder
burden_append<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    mafsummary<-NULL
    mafburden<-NULL
    mafburden_mb<-NULL
    burden<-NULL
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    mafsummary<-as.data.frame(getSampleSummary(maf))
    mafsummary$total_mb<-mafsummary$total/35.8
    mafburden<-median(mafsummary$total)
    mafburden_mb<-median(mafsummary$total_mb)
    burden<-as.data.frame(cbind(mafburden,mafburden_mb))
    burden$Caller<-caller[k]
    burden$Cancer<-cancer[i]
    burden_append<-rbind(burden_append,burden)
  }
}
View(burden_append)

## Merge cancer driver genes with TMB
genes_aggregate$Cohort<-paste(genes_aggregate$Cancer,genes_aggregate$Caller, sep = "_")
burden_append$Cohort<-paste(burden_append$Cancer,burden_append$Caller,sep = "_")
genes_burden<-merge(genes_aggregate,burden_append[,c(1,2,5)],by="Cohort")
dim(genes_burden)
View(genes_burden)
saveRDS(genes_burden,paste(dataDir,"Resources/genes_burden.rds",sep = "/"))

################################################################################
### Step 3: Define reference set (truthset) of cancer driver genes 
################################################################################

## We retrieve our intOGen reference set (truthset) of cancer driver genes (https://www.intogen.org/download)
intogen<-as.data.frame(data.table::fread(paste(dataDir,"Resources/Compendium_Cancer_Genes.tsv",sep = "/")))
head(intogen)
# We subset the driver genes and cancer types (only TCGA cohorts)
intogen<-intogen[,c(1,3,4,10:12)]
intogen<-intogen[grepl("TCGA",intogen$COHORT),]
# Check Cancer Type nomenclature
table(intogen$CANCER_TYPE)
unique(genes_dcast$Cancer[which(!genes_dcast$Cancer %in% intogen$CANCER_TYPE)])
# Modify intOGen Cancer Type nomenclature to match TCGA´s
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="CH"]<-"CHOL"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="MGCT"]<-"TGCT"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="HC"]<-"LIHC"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="RCH"]<-"KICH"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="RCCC"]<-"KIRC"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="DLBCL"]<-"DLBC"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="AML"]<-"LAML"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="RPC"]<-"KIRP"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="ST"]<-"STAD"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="S"]<-"SARC"
intogen$CANCER_TYPE[intogen$CANCER_TYPE=="CM"]<-"SKCM"
# We merge COAD (colon) & READ (rectum) TCGA projects
genes_dcast$Cancer2<-genes_dcast$Cancer
genes_dcast$Cancer2[genes_dcast$Cancer2 == "COAD" | genes_dcast$Cancer2 == "READ"]<-"COREAD"
unique(genes_dcast$Cancer2[which(!genes_dcast$Cancer2 %in% intogen$CANCER_TYPE)])
# We save data
genes_dcast$SYMBOL_CANCER_TYPE<-paste(genes_dcast$SYMBOL,genes_dcast$Cancer,sep = "_")
saveRDS(genes_dcast,paste(dataDir,"Resources/genes_dcast.rds",sep = "/"))
# We set key
intogen$SYMBOL_CANCER_TYPE<-paste(intogen$SYMBOL,intogen$CANCER_TYPE,sep = "_")
# We set a value
intogen$intOGen<-1
# We erase duplicates
intogen<-intogen[!duplicated(intogen$SYMBOL_CANCER_TYPE),]
intogen<-intogen[,-2]
# Duplicate COREAD rows and substitute by COAD and READ respectively
head(intogen)
intogen<-as.data.table(intogen)
intogen<-rbind(intogen,intogen[CANCER_TYPE=="COREAD"][,colnames(intogen)[2] := c("READ")])
intogen$CANCER_TYPE[intogen$CANCER_TYPE == "COREAD"]<-"COAD"
intogen<-as.data.frame(intogen)
intogen$SYMBOL_CANCER_TYPE<-paste(intogen$SYMBOL,intogen$CANCER_TYPE,sep = "_")
head(intogen)
saveRDS(intogen,paste(dataDir,"Resources/intogen_TCGA.rds",sep = "/"))

## We also take the mc3 driver genes list as reference set (truthset) from Bailey et al.,2018
mc3<-read.table(paste(dataDir,'Resources/tissue_driver_pairs.txt',sep = "/"), sep='\t', header=T, fill=T)
colnames(mc3)[1:2]<-c("SYMBOL","CANCER_TYPE")
unique(genes_dcast$Cancer2[which(!genes_dcast$Cancer2 %in% mc3$CANCER_TYPE)])
mc3$CANCER_TYPE[mc3$CANCER_TYPE=="COADREAD"]<-"COREAD"
mc3$SYMBOL_CANCER_TYPE<-paste(mc3$SYMBOL,mc3$CANCER_TYPE,sep = "_")
mc3$mc3<-1
mc3<-mc3[,c(1,2,4:5,13:14)]
mc3<-mc3[!duplicated(mc3$SYMBOL_CANCER_TYPE),]
head(mc3)
#Duplicate COREAD rows and substitute by COAD and READ respectively
head(mc3)
mc3<-as.data.table(mc3)
mc3<-rbind(mc3,mc3[CANCER_TYPE=="COREAD"][,colnames(mc3)[2] := c("READ")])
mc3$CANCER_TYPE[mc3$CANCER_TYPE == "COREAD"]<-"COAD"
mc3<-as.data.frame(mc3)
mc3$SYMBOL_CANCER_TYPE<-paste(mc3$SYMBOL,mc3$CANCER_TYPE,sep = "_")
head(mc3)
#Remove PANCAN
mc3<-mc3[mc3$CANCER_TYPE != "PANCAN",] #235 genes/299 genes with PANCANCER
saveRDS(mc3,paste(dataDir,"Resources/mc3.rds",sep = "/"))

################################################################################
### Step 4: Performance metrics of variant calling strategies to detect cancer driver genes  
################################################################################

# Read data
genes_dcast<-readRDS(paste(dataDir,"Resources/genes_dcast.rds",sep = "/"))
intogen<-readRDS(paste(dataDir,"Resources/intogen_TCGA.rds",sep = "/"))
mc3<-readRDS(paste(dataDir,"Resources/mc3.rds",sep = "/"))
head(genes_dcast)
head(intogen)
head(mc3)
# Merge data at gene and cancer level
metrics<-Reduce(function(x,y) merge(x,y,by="SYMBOL_CANCER_TYPE",all=TRUE) ,list(genes_dcast[,c(3:9,12)],intogen[,-c(1:2)],mc3[,-c(1:2)]))
metrics[is.na(metrics)] <- 0
metrics<-metrics[,c(1:8,12,15,9:11,13:14)]
metrics$Total<-apply(metrics[,2:8],1,function(x) sum(x))
metrics<-metrics[,c(1:8,16,9:15)]
metrics$SYMBOL<-sapply(strsplit(metrics[,1], "\\_"), `[`, 1)
metrics$CANCER<-sapply(strsplit(metrics[,1], "\\_"), `[`, 2)
metrics<-metrics[,c(17:18,1:16)]
head(metrics)
# Save results
saveRDS(metrics,paste(dataDir,"Resources/metrics_TCGA.rds",sep = "/"))
# Write report
subDir<-"Supplementary_files"
dir.create(file.path(dataDir, subDir))
writexl::write_xlsx(metrics,paste(dataDir,"Supplementary_files/Supplementary_file_1_intOGen_results.xlsx",sep = "/"))

## Compute metrics per cancer type

confusion_matrix_metrics_append<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    variant_caller<-NULL
    cancer_type<-NULL
    confusion_matrix_intogen<-NULL
    confusion_matrix_mc3<-NULL
    confusion_matrix_metrics<-NULL
    variant_caller<-caller[k]
    cancer_type<-cancer[i]
    confusion_matrix_intogen<-confusionMatrix(data=factor(metrics[,variant_caller][metrics$CANCER==cancer_type]),reference=factor(metrics$intOGen[metrics$CANCER==cancer_type]), mode = "prec_recall", positive = "1")
    confusion_matrix_mc3<-confusionMatrix(data=factor(metrics[,variant_caller][metrics$CANCER==cancer_type]),reference=factor(metrics$mc3[metrics$CANCER==cancer_type]), mode = "prec_recall", positive = "1")
    confusion_matrix_metrics<-as.data.frame(rbind(confusion_matrix_intogen$byClass[5:7],confusion_matrix_mc3$byClass[5:7]))
    confusion_matrix_metrics$Caller<-caller[k]
    confusion_matrix_metrics$Cancer<-cancer[i]
    confusion_matrix_metrics$Truthset<-c("intOGen","MC3")
    confusion_matrix_metrics_append<-rbind(confusion_matrix_metrics_append,confusion_matrix_metrics)
  }
}

################################################################################
### Step 5: Plot metrics (Figure 2C)  
################################################################################

## Plot performance metrics when detecting cancer driver genes with intOGen
## Plot metrics per cancer type
# Prepare data
confusion_matrix_metrics_append$Caller<-gsub("CONSENSUS2","Consensus2",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("CONSENSUS3","Consensus3",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("MUSE","MuSE",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("MUTECT2","MuTect2",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("SOMATICSNIPER","SomaticSniper",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("UNION","Union",confusion_matrix_metrics_append$Caller)
confusion_matrix_metrics_append$Caller<-gsub("VARSCAN2","VarScan2",confusion_matrix_metrics_append$Caller)
# Subset intOGen and mc3 reference sets
confusion_matrix_metrics_append_intogen<-confusion_matrix_metrics_append[confusion_matrix_metrics_append$Truthset=="intOGen",]
confusion_matrix_metrics_append_mc3<-confusion_matrix_metrics_append[confusion_matrix_metrics_append$Truthset=="MC3",]
head(confusion_matrix_metrics_append_intogen)
head(confusion_matrix_metrics_append_mc3)
summary(confusion_matrix_metrics_append_intogen)
summary(confusion_matrix_metrics_append_mc3)

## Plot metrics per cancer type as boxplots
# intOGen reference set
metrics_boxplot_intogen <- reshape2::melt(confusion_matrix_metrics_append_intogen, id.vars = c("Caller","Truthset","Cancer"),variable.name = "Metrics", value.name = "Score")
metrics_boxplot_intogen$Metrics<-gsub("F1","F1-Score",metrics_boxplot_intogen$Metrics)
metrics_boxplot_intogen$Metrics<-factor(metrics_boxplot_intogen$Metrics, levels = c("Precision","Recall","F1-Score"))
head(metrics_boxplot_intogen)
# Plot Figure 2C (left panel)
m_boxplot_intogen<-ggplot(data = metrics_boxplot_intogen, aes(reorder_within(Caller,Score,Metrics,mean),Score,)) + theme_bw(base_size = 12) + geom_boxplot(aes(fill = Caller),width = 0.8, outlier.colour=NA) + scale_x_reordered() + ylab ("Score") + xlab (element_blank()) + facet_grid (as.factor(metrics_boxplot_intogen$Truthset)~as.factor (metrics_boxplot_intogen$Metrics), scales = "free_x") + theme(axis.text.x = element_text (angle = 45, hjust = 1)) + theme(aspect.ratio = 1, legend.position = "bottom", legend.margin=margin(t=-15),axis.ticks.x=element_blank(), axis.text.x = element_blank(),legend.title = element_blank(),legend.key.size = unit(1,"line"), legend.text = element_text(size = 12, colour = "black", face = "bold")) + guides(fill = guide_legend(nrow = 2)) + scale_fill_manual(values=c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) + scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + theme(strip.background = element_rect(colour="black", fill="antiquewhite",size=1, linetype="solid")) + theme(axis.text.y = element_text(color="black", size=12, face="bold"),axis.title.y = element_text(color="black", size=14, face="bold"),strip.text.x = element_text(size=12, color="black",face="bold"),strip.text.y = element_text(size=12, color="black",face="bold")) 
ggsave(plot = m_boxplot_intogen, width = 7, height = 7, dpi = 300,filename = paste(dataDir,"Figures/Figure2C_intogen.png",sep = "/"))
saveRDS(m_boxplot_intogen,paste(dataDir,"Figures/Figure2C_intogen.rds",sep = "/"))
# mc3 reference set
metrics_boxplot_mc3 <- reshape2::melt(confusion_matrix_metrics_append_mc3, id.vars = c("Caller","Truthset","Cancer"),variable.name = "Metrics", value.name = "Score")
metrics_boxplot_mc3$Metrics<-gsub("F1","F1-Score",metrics_boxplot_mc3$Metrics)
metrics_boxplot_mc3$Metrics<-factor(metrics_boxplot_mc3$Metrics, levels = c("Precision","Recall","F1-Score"))
head(metrics_boxplot_mc3)
# Plot Figure 2C (right panel)
m_boxplot_mc3<-ggplot(data = metrics_boxplot_mc3, aes(reorder_within(Caller,Score,Metrics,mean),Score,)) + theme_bw(base_size = 12) + geom_boxplot(aes(fill = Caller),width = 0.8, outlier.colour=NA) + scale_x_reordered() + ylab ("Score") + xlab (element_blank()) + facet_grid (as.factor(metrics_boxplot_mc3$Truthset)~as.factor (metrics_boxplot_mc3$Metrics), scales = "free_x") + theme(axis.text.x = element_text (angle = 45, hjust = 1)) + theme(aspect.ratio = 1, legend.position = "bottom", legend.margin=margin(t=-15),axis.ticks.x=element_blank(), axis.text.x = element_blank(),legend.title = element_blank(),legend.key.size = unit(1,"line"), legend.text = element_text(size = 12, colour = "black", face = "bold")) + guides(fill = guide_legend(nrow = 2)) + scale_fill_manual(values=c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) + scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1)) + theme(strip.background = element_rect(colour="black", fill="antiquewhite",size=1, linetype="solid")) + theme(axis.text.y = element_text(color="black", size=12, face="bold"),axis.title.y = element_text(color="black", size=14, face="bold"),strip.text.x = element_text(size=12, color="black",face="bold"),strip.text.y = element_text(size=12, color="black",face="bold")) 
ggsave(plot = m_boxplot_mc3, width = 7, height = 7, dpi = 300,filename = paste(dataDir,"Figures/Figure2C_mc3.png",sep = "/"))
saveRDS(m_boxplot_mc3,paste(dataDir,"Figures/Figure2C_mc3.rds",sep = "/"))
# Write report
metrics_boxplot_intogen_mc3_per_cancer<-rbind(metrics_boxplot_intogen,metrics_boxplot_mc3)
writexl::write_xlsx(metrics_boxplot_intogen_mc3_per_cancer,paste(dataDir,"Supplementary_files/Supplementary_file_2_metrics_boxplot_intogen_mc3_per_cancer.xlsx",sep = "/"))

################################################################################
### Step 6: Plot best performing variant calling strategy (Figure 2D)  
################################################################################

# Prepare data for alluvial plot
truthset<-unique(confusion_matrix_metrics_append$Truthset)
alluvial<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  alluvial_cohort<-NULL
  alluvial_cohort<-confusion_matrix_metrics_append[confusion_matrix_metrics_append$Cancer==cancer[i],]
  for(z in 1:length(truthset)){
    print(paste(i,cancer[i],truthset[z],sep = "_"))
    alluvial_cohort_truthset<-NULL
    alluvial_cohort_truthset_final<-NULL
    alluvial_cohort_truthset<-alluvial_cohort[alluvial_cohort$Truthset == truthset[z],]
    alluvial_cohort_truthset_final<-as.data.frame(cbind(alluvial_cohort_truthset$Caller[which.max(alluvial_cohort_truthset$Precision)],alluvial_cohort_truthset$Caller[which.max(alluvial_cohort_truthset$Recall)],alluvial_cohort_truthset$Caller[which.max(alluvial_cohort_truthset$F1)],truthset[z],cancer[i]))
    colnames(alluvial_cohort_truthset_final)<-c("Precision","Recall","F1-Score","Truthset","Cancer")
    alluvial<-rbind(alluvial,alluvial_cohort_truthset_final)
  }
}
#Save data
saveRDS(alluvial,paste(dataDir,"Resources/alluvial.rds",sep = "/"))
head(alluvial)

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

## Alluvial plot (Figure 2D)
# Read data
alluvial<-readRDS(paste(dataDir,"Resources/alluvial.rds",sep = "/"))
# Prepare data
alluvial$Group<-ifelse(alluvial$Cancer %in% Hematologic_lymphatic,"Hematologic and lymphatic",ifelse(alluvial$Cancer %in% Gynecologic, "Gynecologic", ifelse(alluvial$Cancer %in% Urologic, "Urologic", ifelse(alluvial$Cancer %in% Endocrine, "Endocrine",ifelse(alluvial$Cancer %in% Gastrointestinal,"Gastrointestinal",ifelse(alluvial$Cancer %in% Thoracic,"Thoracic",ifelse(alluvial$Cancer %in% CNS,"Central Nervous System",ifelse(alluvial$Cancer %in% Soft_tissue,"Soft tissue","Other"))))))))
alluvial$Group<-factor(alluvial$Group,levels = c("Central Nervous System","Endocrine","Gastrointestinal","Gynecologic","Hematologic and lymphatic","Soft tissue","Thoracic","Urologic","Other"))
alluvial_intogen<-alluvial[alluvial$Truthset=="intOGen",]
alluvial_mc3<-alluvial[alluvial$Truthset=="MC3",]
# mc3 reference set
alluvial_mc3_plot<-ggplot(data = alluvial_mc3,
                      aes(axis2 = Group, axis1 = `F1-Score`)) +
                      scale_x_discrete(limits = c("F1-Score","Group"), expand = c(.0, .0)) +
                      geom_alluvium(aes(fill = `F1-Score`)) +
                      geom_stratum(color = "black") +
                      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                      theme_classic() +
                      scale_fill_manual(values = c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                      ggtitle("Best Caller")
ggsave(plot = alluvial_mc3_plot, width = 15.4, height = 5.5, dpi = 300, filename = paste(dataDir,"Figures/Figure2D_mc3.png",sep = "/"))
# intOGen reference set
alluvial_intogen_plot<-ggplot(data = alluvial_intogen,
                              aes(axis2 = Group, axis1 = `F1-Score`)) +
                              scale_x_discrete(limits = c("F1-Score","Group"), expand = c(.0, .0)) +
                              geom_alluvium(aes(fill = `F1-Score`)) +
                              geom_stratum(color = "black") +
                              geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                              theme_classic() +
                              scale_fill_manual(values = c("#8F7700FF","#0073C2FF","#EFC000FF","#868686FF","#CD534CFF","#7AA6DCFF","#003C67FF")) +
                              ggtitle("Best Caller intOGen")
ggsave(plot = alluvial_intogen_plot, width = 9.5, height = 5, dpi = 300, filename = paste(dataDir,"Figures/Figure2D_intogen.png",sep = "/"))
# Save results
saveRDS(alluvial_intogen_plot,paste(dataDir,"Figures/Figure2D_intogen.rds",sep = "/"))
saveRDS(alluvial_mc3_plot,paste(dataDir,"Figures/Figure2D_mc3.rds",sep = "/"))


################################################################################
### Step 7: Mutations in cancer driver genes 
################################################################################

# Read data
metrics<-readRDS(paste(dataDir,"metrics_TCGA.rds",sep = "/"))
head(metrics)
# Select cancer driver genes in intOGen & mc3
table(metrics$intOGen,metrics$mc3)
metrics_drivers_intogen_mc3<-metrics[metrics$intOGen == 1 & metrics$mc3 == 1,]
dim(metrics_drivers_intogen_mc3)
# Select cancer driver genes in more than one tumor type
drivers_duplicated<-unique(metrics_drivers_intogen_mc3$SYMBOL[duplicated(metrics_drivers_intogen_mc3$SYMBOL)])
metrics_drivers_intogen_mc3_duplicates<-metrics_drivers_intogen_mc3[metrics_drivers_intogen_mc3$SYMBOL %in% drivers_duplicated,]
# Select cancer driver genes called by at least one method
head(metrics_drivers_intogen_mc3_duplicates)
table(metrics_drivers_intogen_mc3_duplicates$Total)
metrics_drivers_intogen_mc3_duplicates_method<-metrics_drivers_intogen_mc3_duplicates[metrics_drivers_intogen_mc3_duplicates$Total > 0,]
dim(metrics_drivers_intogen_mc3_duplicates_method)
# Check cancer driver gene roles
table(metrics_drivers_intogen_mc3_duplicates_method$ROLE,metrics_drivers_intogen_mc3_duplicates_method$Tumor.suppressor.or.oncogene.prediction..by.20.20..)
# Check final cancer driver genes list
driver_genes_somatic_mutations<-unique(metrics_drivers_intogen_mc3_duplicates_method$SYMBOL)
head(metrics_drivers_intogen_mc3_duplicates_method)
table(metrics_drivers_intogen_mc3_duplicates_method$CANCER)
table(metrics_drivers_intogen_mc3_duplicates_method$SYMBOL)
# Select all cancer driver genes in intogen & mc3
all_drivers_intogen_mc3<-unique(metrics_drivers_intogen_mc3$SYMBOL)

## Retrieve all somatic mutations in cancer driver genes
driver_mutations<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    mafsummary<-NULL
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    mafsummary<-as.data.frame(getGeneSummary(maf))
    mafsummary<-mafsummary[mafsummary$Hugo_Symbol %in% all_drivers_intogen_mc3,]
    mafsummary<-mafsummary[,c("Hugo_Symbol","Missense_Mutation","Nonsense_Mutation","MutatedSamples")]
    mafsummary$total_miss_non_sense_mutations<-mafsummary$Missense_Mutation + mafsummary$Nonsense_Mutation
    mafsummary$Caller<-caller[k]
    mafsummary$Cancer<-cancer[i]
    driver_mutations<-rbind(driver_mutations,mafsummary)
  }
}
head(driver_mutations)

## Retrieve cohort size
cohort_size<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  maf<-NULL
  mafsummary<-NULL
  maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_","UNION",".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
  maf<-as.data.frame(maf@summary)
  mafsummary<-maf[3,1:2]
  mafsummary$Cancer<-cancer[i]
  mafsummary<-mafsummary[,2:3]
  colnames(mafsummary)[1]<-"Samples"
  rownames(mafsummary)<-NULL
  cohort_size<-rbind(cohort_size,mafsummary)
}
colnames(cohort_size)<-c("Samples","Cancer")
cohort_size$Samples<-as.numeric(cohort_size$Samples)
head(cohort_size)

## Create data frame with somatic mutations information
driver_mutations_size<-merge(driver_mutations,cohort_size, by = "Cancer")
driver_mutations_size$Percentage<-round((driver_mutations_size$MutatedSamples/driver_mutations_size$Samples)*100,2)
driver_mutations_size$SYMBOL_CANCER_TYPE<-paste(driver_mutations_size$Hugo_Symbol,driver_mutations_size$Cancer,sep = "_")
driver_mutations_size$Percentage<-driver_mutations_size$Percentage/100
head(metrics_drivers_intogen_mc3_duplicates_method)
head(driver_mutations_size)
# Select only driver genes in each cancer type
driver_mutations_size_cohort<-driver_mutations_size[driver_mutations_size$SYMBOL_CANCER_TYPE %in% metrics_drivers_intogen_mc3_duplicates_method$SYMBOL_CANCER_TYPE,]
head(driver_mutations_size_cohort)
dim(driver_mutations_size_cohort)
table(driver_mutations_size_cohort$Hugo_Symbol)
table(driver_mutations_size_cohort$Cancer)

## Plot mutations in cancer driver genes
#Prepare data
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "CONSENSUS2"]<-"Consensus2"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "CONSENSUS3"]<-"Consensus3"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "MUSE"]<-"MuSE"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "MUTECT2"]<-"MuTect2"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "SOMATICSNIPER"]<-"SomaticSniper"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "UNION"]<-"Union"
driver_mutations_size_cohort$Caller[driver_mutations_size_cohort$Caller == "VARSCAN2"]<-"VarScan2"
driver_mutations_size_cohort$Cancer_size<-paste0(driver_mutations_size_cohort$Cancer," (n=",driver_mutations_size_cohort$Samples,")")
# Save data
saveRDS(driver_mutations_size_cohort,paste(dataDir,"Resources/driver_mutations_size_cohort.rds",sep = "/"))

################################################################################
### Step 8: Plot Supplementary Figure S5
################################################################################

# Read data
driver_mutations_size_cohort<-readRDS(paste(dataDir,"Resources/driver_mutations_size_cohort.rds",sep = "/"))
# Select most frequently mutated cancer driver genes
genes_interest<-c("TP53","KRAS","PTEN","PIK3CA")
driver_mutations_size_genes<-driver_mutations_size[driver_mutations_size$Hugo_Symbol %in% genes_interest,]
table(driver_mutations_size_genes$Hugo_Symbol)
saveRDS(driver_mutations_size_genes,paste(dataDir,"Figures/FigureS5.rds",sep = "/"))
# Plot Figure S5
driver_mutations_size_genes<-readRDS(paste(dataDir,"Figures/FigureS5.rds",sep = "/"))
for(i in 1:length(genes_interest)){
  print(i)
  print(genes_interest[i])
  plot_mutations<-NULL
  plot_mutations<-ggplot(data = driver_mutations_size_genes[driver_mutations_size_genes$Hugo_Symbol == genes_interest[i],], aes(x= Cancer, y = Caller, fill=Percentage)) + scale_fill_viridis(na.value = "white",option="A", begin = 0, end = 1, direction = -1, label = scales::percent_format(accuracy = 1)) + geom_tile(color="gray20",size=.7, stat = "identity") + facet_wrap (~as.factor (driver_mutations_size_genes[driver_mutations_size_genes$Hugo_Symbol == genes_interest[i],]$Hugo_Symbol)) + theme(axis.ticks=element_blank()) + theme (axis.text.x = element_text (angle = 90, vjust=.5, hjust=1, size = 10, colour = "black"),axis.text.y = element_text(size = 10, colour = "black")) + xlab(element_blank()) + ylab(element_blank()) + ggtitle(element_blank()) + labs(fill="Samples") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), strip.background = element_rect(colour="black", fill="antiquewhite",size=1, linetype="solid"),legend.position = "top",legend.key.size = unit(1,"line"), legend.key.height =  unit(.75,"line"), legend.key.width =  unit(1.25,"line"), legend.text = element_text(size = 9, colour = "black"), legend.title = element_text(size = 9, color = "black"), strip.text.x = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9),legend.justification = "center",legend.margin = margin(0.0,0.0),legend.box.margin = margin(-10,-10,-10,-10))
  ggsave(plot = plot_mutations, width = 5, height = 2.5, dpi = 300, filename = paste0(dataDir,"/Figures/",genes_interest[i],"_","FigureS5.png"))
}

################################################################################
### Step 9: Plot Figures 2A and 2B
################################################################################

# Read data
genes_burden<-readRDS(paste(dataDir,"Resources/genes_burden.rds",sep = "/"))
driver_mutations_size_cohort<-readRDS(paste(dataDir,"Resources/driver_mutations_size_cohort.rds",sep = "/"))
head(genes_burden)
head(driver_mutations_size_cohort)
# Add cohort size to genes burden 
genes_burden_size<-merge(genes_burden,unique(driver_mutations_size_cohort[,c(1,8,11)]), by = "Cancer")
head(genes_burden_size)
saveRDS(genes_burden_size,paste(dataDir,"Resources/genes_burden_size.rds",sep = "/"))

## Plot Figures 2A and 2B
genes_burden<-readRDS(paste(dataDir,"Resources/genes_burden_size.rds",sep = "/"))
pancancer_colors<-c("#67000D","#A50F15","#EF3B2C","#FC9272","#FEE0D2","#BCBDDC","#807DBA",
                    "#54278F","#3F007D","#08306B","#08519C","#4292C6","#9ECAE1","#000000",
                    "#525252","#969696","#BDBDBD","#D9D9D9","#80CDC1","#35978F","#01665E",
                    "#006D2C","#41AB5D","#A1D99B","#FEFF4F","#FED976","#FD8D3C","#8C510A",
                    "#BF812D","#DFC27D","#FA9FB5","#F768A1","#DD3497")
#Figure2A
Figure2<-ggplot(data = genes_burden, aes(y = Key, x = mafburden_mb)) + theme_bw(base_size = 12) + geom_point(aes(colour = Cancer), size = 2)  + ylab ("intOGen detected driver genes") + xlab ("log10 (median number of mutations per megabase)") + stat_cor(method = "spearman", alternative = "two.sided", cor.coef.name = "rho", label.y = 88, size = 5) + stat_smooth(method = "lm", colour = "black", size=0.5, se = T) + scale_x_log10() + scale_y_continuous(breaks=seq(0,90,10), limits = c(0,90)) + scale_color_manual(values=pancancer_colors) + theme(aspect.ratio = 1, legend.position = "right", legend.background = element_blank(), legend.title = element_blank()) + theme(axis.title = element_text(size = 14, face = "bold"),axis.text = element_text(size = 12))
ggsave(plot = Figure2, width = 6.5, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/Figure2A.png",sep = "/"))
# Save data
saveRDS(Figure2,paste(dataDir,"Figures/Figure2A.rds",sep = "/"))
# Figure2B
set.seed(24)
Figure2_SampleSize<-ggplot(data = genes_burden, aes(y = Key, x = Samples)) + theme_bw(base_size = 12) + geom_point(aes(colour = Cancer), size = 2, position=position_jitter(h=0, w=0.01))  + ylab ("intOGen detected driver genes") + xlab ("Number of samples") + stat_cor(method = "spearman", alternative = "two.sided", cor.coef.name = "rho", label.y = 88, size = 5) + stat_smooth(method = "lm", colour = "black", size=0.5, se = T) + scale_y_continuous(breaks=seq(0,90,10), limits = c(0,90)) + scale_color_manual(values=pancancer_colors) + theme(aspect.ratio = 1, legend.position = "right", legend.background = element_blank(), legend.title = element_blank()) + theme(axis.title = element_text(size = 14, face = "bold"),axis.text = element_text(size = 12)) + scale_x_log10()
ggsave(plot = Figure2_SampleSize, width = 6.5, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/Figure2B.png",sep = "/"))
# Save data
saveRDS(Figure2_SampleSize,paste(dataDir,"Figures/Figure2B.rds",sep = "/"))
#===============================================================================