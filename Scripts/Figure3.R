#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure 3
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
### Step 1: Retrieve all somatic mutations in all cancer driver genes
################################################################################

# Read performance metrics 
metrics<-readRDS(paste(dataDir,"Resources/metrics_TCGA.rds",sep = "/"))
head(metrics)
# Select all cancer driver genes (defined by intOGen & mc3)
table(metrics$intOGen,metrics$mc3)
metrics_drivers_intogen_or_mc3<-metrics[metrics$intOGen == 1 | metrics$mc3 == 1,]
dim(metrics_drivers_intogen_or_mc3)
unique(metrics_drivers_intogen_or_mc3$SYMBOL)
head(metrics_drivers_intogen_or_mc3)
# Read all somatic mutations data
maf_append_purity_ready_clinical<-readRDS(paste(dataDir,"Resources/maf_append_purity_ready_clinical.rds",sep = "/"))
head(maf_append_purity_ready_clinical)
dim(maf_append_purity_ready_clinical)
length(unique(maf_append_purity_ready_clinical$name))
table(maf_append_purity_ready_clinical$Cancer)
# Select somatic mutations in all cancer driver genes
maf_append_purity_ready_clinical_drivers<-maf_append_purity_ready_clinical[maf_append_purity_ready_clinical$Hugo_Symbol %in% unique(metrics_drivers_intogen_or_mc3$SYMBOL),]
dim(maf_append_purity_ready_clinical_drivers)
unique(metrics_drivers_intogen_or_mc3$SYMBOL)[which(!unique(metrics_drivers_intogen_or_mc3$SYMBOL) %in% unique(maf_append_purity_ready_clinical_drivers$Hugo_Symbol))]
# Create Symbol_Cancer_type key column
maf_append_purity_ready_clinical_drivers$SYMBOL_CANCER_TYPE<-paste(maf_append_purity_ready_clinical_drivers$Hugo_Symbol,maf_append_purity_ready_clinical_drivers$Cancer,sep = "_")
head(maf_append_purity_ready_clinical_drivers)
# Select only those cancer driver genes somatic mutations found in the same cancer types where their respective cancer driver genes have been defined
maf_append_purity_ready_clinical_drivers_cancer<-maf_append_purity_ready_clinical_drivers[maf_append_purity_ready_clinical_drivers$SYMBOL_CANCER_TYPE %in% unique(metrics_drivers_intogen_or_mc3$SYMBOL_CANCER_TYPE),]
dim(maf_append_purity_ready_clinical_drivers_cancer)

################################################################################
### Step 2: Prepare data for plot
################################################################################

# Put caller column as last column
table(maf_append_purity_ready_clinical_drivers_cancer$Cancer,maf_append_purity_ready_clinical_drivers_cancer$type)
colnames(maf_append_purity_ready_clinical_drivers_cancer)
maf_append_purity_ready_clinical_drivers_cancer<-maf_append_purity_ready_clinical_drivers_cancer[,c(1:134,137:170,135:136)]
# Select missense and nonsense mutations
table(maf_append_purity_ready_clinical_drivers_cancer$Variant_Classification)
maf_append_purity_ready_clinical_drivers_cancer<-maf_append_purity_ready_clinical_drivers_cancer[maf_append_purity_ready_clinical_drivers_cancer$Variant_Classification=="Missense_Mutation" | maf_append_purity_ready_clinical_drivers_cancer$Variant_Classification=="Nonsense_Mutation",]
dim(maf_append_purity_ready_clinical_drivers_cancer)
#Transform blank spaces to NA
maf_append_purity_ready_clinical_drivers_cancer[maf_append_purity_ready_clinical_drivers_cancer==""] <- NA

# Reshape data
maf_append_purity_ready_clinical_drivers_cancer_dcast<-reshape2::dcast(maf_append_purity_ready_clinical_drivers_cancer, name+MC3_Overlap+t_vaf_purity_ploidy_final_mean+total_coverage_mean+Variant_Classification+SYMBOL_CANCER_TYPE~Caller, fun.aggregate = length)
head(maf_append_purity_ready_clinical_drivers_cancer_dcast)
colnames(maf_append_purity_ready_clinical_drivers_cancer_dcast)[c(1:13)]<-c("Variant","MC3_Overlap","VAF","Coverage","Variant Classification","SYMBOL_CANCER_TYPE","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
table(maf_append_purity_ready_clinical_drivers_cancer_dcast$MC3) #95.86% overlap MC3

## Merge cancer driver genes data
# Select information in maf_append_purity_ready_clinical_drivers_cancer_dcast file
metrics_drivers_intogen_or_mc3_subset<-metrics_drivers_intogen_or_mc3[metrics_drivers_intogen_or_mc3$SYMBOL_CANCER_TYPE %in% maf_append_purity_ready_clinical_drivers_cancer_dcast$SYMBOL_CANCER_TYPE,]
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged<-merge(maf_append_purity_ready_clinical_drivers_cancer_dcast,metrics_drivers_intogen_or_mc3_subset,by="SYMBOL_CANCER_TYPE",all = TRUE)
dim(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged)
dim(maf_append_purity_ready_clinical_drivers_cancer_dcast)
head(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged)
# Prepare for upset plot
callers_drivers<-c("MC3_Overlap",callers)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$MC3_Overlap<-ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$MC3_Overlap == TRUE, 1, 0)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged[c(callers_drivers)]<-maf_append_purity_ready_clinical_drivers_cancer_dcast_merged[c(callers_drivers)] != 0
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$`Variant Classification`<-gsub("Missense_Mutation","Missense",maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$`Variant Classification`)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$`Variant Classification`<-gsub("Nonsense_Mutation","Nonsense",maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$`Variant Classification`)
head(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged)
# Set oncogene or tsg function
table(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$ROLE,maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Role<-ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..=="oncogene" | maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..=="possible oncogene","Oncogene",ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..=="tsg" | maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..=="possible tsg","Tumor suppressor gene",ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$ROLE=="Act","Oncogene",ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$ROLE=="LoF","Tumor suppressor gene","Unknown"))))
head(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged)
table(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$ROLE,maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Tumor.suppressor.or.oncogene.prediction..by.20.20..,maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Role)
table(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Role)
# Check coverage information
summary(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage_level<-ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage<45,"<45",
                                                                                    ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage<=75,"45-75",
                                                                                           ifelse(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage<=175,">75-135",">135")))
head(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged)
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage_level<-factor(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage_level,levels = c(">135",">75-135","45-75","<45"))
table(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged$Coverage_level)
#Save data
saveRDS(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged,paste(dataDir,"Resources/maf_append_purity_ready_clinical_drivers_cancer_dcast_merged.rds",sep = "/"))

################################################################################
### Step 3: Plot Figure 3
################################################################################

# Read data
maf_append_purity_ready_clinical_drivers_cancer_dcast_merged<-readRDS(paste(dataDir,"Resources/maf_append_purity_ready_clinical_drivers_cancer_dcast_merged.rds",sep = "/"))
# Set color palette
cbPalette <- c("#E69F00","#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Upset plot
Figure_driver_mutations<-
  (
    upset(maf_append_purity_ready_clinical_drivers_cancer_dcast_merged, callers, name=NULL, width_ratio=0.2, height_ratio = 0.8,  intersections = list(c("Union","MuTect2","VarScan2","MuSE","SomaticSniper","Consensus2","Consensus3"),c("Union","MuTect2","VarScan2","MuSE","Consensus2","Consensus3"),c("Union","MuTect2"),c("Union","VarScan2","MuSE","SomaticSniper","Consensus2","Consensus3"),c("Union","MuTect2","MuSE","Consensus2"),c("Union","VarScan2"),c("Union","MuTect2","MuSE","SomaticSniper","Consensus2","Consensus3"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuSE"),c("Union","SomaticSniper")), keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("Variant Classification"=(ggplot(mapping=aes(fill=`Variant Classification`))+geom_bar(stat='count', position='fill')+theme(legend.key.size = unit(0.85,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name="Variant",values=c('Missense'="lightgoldenrod", 'Nonsense'="deepskyblue4"))+ylab("Variant Classification")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "VAF"=(ggplot(mapping=aes(y=VAF))+geom_jitter(aes(color=Coverage_level),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(0.85,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>135'="cadetblue1", '>75-135'="aquamarine1", '45-75'="antiquewhite3", '<45'="antiquewhite")))),
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),text=list(size=3),mapping=aes(fill=Role)) + scale_fill_manual(values=c('Oncogene'="darkgoldenrod1", 'Tumor suppressor gene'="darkblue", "Unknown"="darkgoldenrod4")) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Role)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=38000)+theme(legend.position = "none")+scale_fill_manual(values=c("Oncogene"="darkgoldenrod1","Tumor suppressor gene"="darkblue","Unknown"="darkgoldenrod4"))))
    + ggtitle("Mutations in cancer driver genes") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_driver_mutations, width = 9, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/Figure3.png",sep = "/"))
# Save data
saveRDS(Figure_driver_mutations,paste(dataDir,"Figures/Figure3.rds",sep="/"))
#===============================================================================