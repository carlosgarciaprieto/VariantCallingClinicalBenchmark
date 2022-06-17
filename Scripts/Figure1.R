#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure 1
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
### Step 1: data loading
################################################################################

## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          "LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          "THCA","THYM","UCEC","UCS","UVM")

caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")

## Read all mafs
# For each of the 33 TCGA cancer types, we downloaded the four different Somatic aggregated MAF files with all the somatic mutations for each variant caller (MuSE, MuTect2, SomaticSniper and VarScan2), created the Consensus of 2, Consensus of 3 and Union files, and saved them under Maf_files folder
subDir<-"Maf_files"
dir.create(file.path(dataDir, subDir))
maf_append<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    maf@data$name<-paste(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode, sep = ":")
    maf@data$Caller<-caller[k]
    maf@data$Cancer<-cancer[i]
    maf_append<-rbind(maf_append,maf@data)
  }
}
# Explore data
View(maf_append)

################################################################################
### Step 2: data processing (adjust VAF by purity and ploidy)
################################################################################

# Transform blank spaces to NA!!!
maf_append[maf_append==""] <- NA
# Put Caller column as last column
head(maf_append)
colnames(maf_append)
maf_append<-maf_append[,c(1:121,123:122)]
# Compute Variant Allele Frequency (VAF)
maf_append$t_vaf<-maf_append$t_alt_count/(maf_append$t_ref_count+maf_append$t_alt_count)

## Adjust VAF by purity and ploidy
# Read ABSOLUTE purity and ploidy data downloaded from MC3 (https://gdc.cancer.gov/about-data/publications/pancanatlas) into Resources folder
subDir2<-"Resources"
dir.create(file.path(dataDir, subDir2))
purity<-data.table::fread(paste(dataDir,"Resources/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",sep = "/"))
purity$name<-gsub('.{5}$', '', purity$sample)
purity$name<-substr(purity$sample,1,20)
head(purity)

## Merge maf file with purity file
# Merge by exact TCGA barcode sample
length(which(unique(maf_append$Tumor_Sample_Barcode) %in% unique(purity$sample))) #Only 19 samples
# Merge by plate TCGA barcode
head(maf_append)
maf_append$PurityID<-substr(maf_append$Tumor_Sample_Barcode,1,20)
# Check overlap
length(which(unique(maf_append$PurityID) %in% unique(purity$name))) #8673 samples (85% samples)
# Check samples with purity/ploidy info
length(which(unique(maf_append$PurityID) %in% unique(purity$array))) #9871 samples with purity/ploidy info
# Check total number of samples
length(unique(maf_append$Tumor_Sample_Barcode)) #10189 total samples
# Select ids in purity file present in maf file
purity_filter<-subset(purity, name %in% maf_append$PurityID)
dim(purity)
dim(purity_filter)
colnames(purity_filter)
# Merge maf with purity/ploidy
maf_append_purity<-merge(maf_append,purity_filter[,c(11,4:6,8:9)],by.x="PurityID",by.y="name",all = TRUE)
dim(maf_append_purity)
dim(maf_append)
head(maf_append_purity)

## Compute VAF adjusted by purity and ploidy
# Adjust purity
maf_append_purity$t_vaf_purity<-maf_append_purity$t_vaf/maf_append_purity$`Cancer DNA fraction`
length(which(maf_append_purity$t_vaf_purity>1)) #341055 variants
length(which(maf_append_purity$t_vaf_purity>1))/(dim(maf_append_purity)[1])*100 #1.71% of variants
length(which(maf_append_purity$t_vaf_purity<1)) #16482437 variants
summary(maf_append_purity$t_vaf_purity)
summary(maf_append_purity$t_vaf)
# Adjust ploidy
table(maf_append_purity$`Genome doublings`)
summary(maf_append_purity$ploidy)
maf_append_purity$t_vaf_purity_ploidy<-maf_append_purity$t_vaf_purity/(maf_append_purity$ploidy/2)
length(which(round(maf_append_purity$t_vaf_purity_ploidy,1)>1)) #176243 variants
length(which(round(maf_append_purity$t_vaf_purity_ploidy,1)>1))/(dim(maf_append_purity)[1])*100 #0.88% of variants
length(which(maf_append_purity$t_vaf_purity_ploidy<1)) #16507666 variants
summary(maf_append_purity$t_vaf_purity_ploidy)
summary(maf_append_purity$t_vaf_purity)
summary(maf_append_purity$t_vaf)
head(maf_append_purity)
# Finally replace VAFs > 1 (samples are matched by plate, no exact matches) and NAs with unadjusted VAFs
sum(is.na(maf_append_purity$t_vaf_purity_ploidy))/(dim(maf_append_purity)[1])*100 #15.6% NAs
maf_append_purity$t_vaf_purity_ploidy_final<-ifelse(is.na(maf_append_purity$t_vaf_purity_ploidy), maf_append_purity$t_vaf,ifelse(round(maf_append_purity$t_vaf_purity_ploidy,1) >1, maf_append_purity$t_vaf,maf_append_purity$t_vaf_purity_ploidy))
dim(maf_append_purity)
sum(is.na(maf_append_purity$t_vaf)) #2 samples missing t_vaf info
head(maf_append_purity)

################################################################################
### Step 3: reshape data and compute mean VAF and coverage by variant
################################################################################

# Compute mean vaf and coverage by variant
t_vaf_purity_ploidy_final_mean<-aggregate(maf_append_purity$t_vaf_purity_ploidy_final, list(maf_append_purity$name), FUN=mean)
maf_total_coverage_mean<-aggregate(maf_append_purity$t_ref_count+maf_append_purity$t_alt_count, list(maf_append_purity$name), FUN=mean)
colnames(t_vaf_purity_ploidy_final_mean)<-c("name","t_vaf_purity_ploidy_final_mean")
colnames(maf_total_coverage_mean)<-c("name","total_coverage_mean")
maf_total_coverage_mean$total_coverage_mean<-round(maf_total_coverage_mean$total_coverage_mean)
head(t_vaf_purity_ploidy_final_mean)
dim(t_vaf_purity_ploidy_final_mean)
length(unique(maf_append_purity$name))
head(maf_total_coverage_mean)
dim(maf_total_coverage_mean)
# Merge both means
t_vaf_purity_ploidy_final_coverage_mean<-merge(t_vaf_purity_ploidy_final_mean,maf_total_coverage_mean,by="name")
head(t_vaf_purity_ploidy_final_coverage_mean)
dim(t_vaf_purity_ploidy_final_coverage_mean)

## Merge with maf file
maf_append_purity_ready<-merge(maf_append_purity,t_vaf_purity_ploidy_final_coverage_mean, by = "name")
dim(maf_append_purity_ready)
dim(maf_append_purity)
# Transform blank spaces to NA!!!
maf_append_purity_ready[maf_append_purity_ready==""] <- NA
head(maf_append_purity_ready)

################################################################################
### Step 4: Add clinical information
################################################################################

## Read clinical information downloaded from MC3 (https://gdc.cancer.gov/about-data/publications/pancanatlas) into Resources folder
clinical<-as.data.frame(readxl::read_excel(paste(dataDir,"Resources/TCGA-CDR-SupplementalTableS1.xlsx",sep = "/"),sheet = 1))[,-1]
head(clinical)
dim(clinical)
table(clinical$ajcc_pathologic_tumor_stage)
table(clinical$clinical_stage)
# Add column with summarize Stage information: Stages I-IV (otherwise considered Not Available)
clinical$ajcc_stage_parsed<-ifelse(clinical$ajcc_pathologic_tumor_stage == "Stage I" | clinical$ajcc_pathologic_tumor_stage == "Stage IA" | clinical$ajcc_pathologic_tumor_stage == "Stage IB", "Stage I",
                                   ifelse(clinical$ajcc_pathologic_tumor_stage == "Stage II" | clinical$ajcc_pathologic_tumor_stage == "Stage IIA" | clinical$ajcc_pathologic_tumor_stage == "Stage IIB" | clinical$ajcc_pathologic_tumor_stage == "Stage IIC", "Stage II",
                                          ifelse(clinical$ajcc_pathologic_tumor_stage == "Stage III" | clinical$ajcc_pathologic_tumor_stage == "Stage IIIA" | clinical$ajcc_pathologic_tumor_stage == "Stage IIIB" | clinical$ajcc_pathologic_tumor_stage == "Stage IIIC", "Stage III",
                                                 ifelse(clinical$ajcc_pathologic_tumor_stage == "Stage IV" | clinical$ajcc_pathologic_tumor_stage == "Stage IVA" | clinical$ajcc_pathologic_tumor_stage == "Stage IVB" | clinical$ajcc_pathologic_tumor_stage == "Stage IVC", "Stage IV",
                                                        "Not Available"))))

table(clinical$ajcc_pathologic_tumor_stage,clinical$ajcc_stage_parsed)

## Merge maf file with clinical information
# Prepare maf files data
maf_append_purity_ready$MC3_Overlap<-toupper(maf_append_purity_ready$MC3_Overlap)
table(maf_append_purity_ready$MC3_Overlap)
# Put caller column as last column
head(maf_append_purity_ready)
colnames(maf_append_purity_ready)
maf_append_purity_ready<-maf_append_purity_ready[,c(1:122,125:135,124:123)]
table(maf_append_purity_ready$GDC_Validation_Status) #703270 Validated variants 3.5% 

## Format clinical data
# Create PatientID column
maf_append_purity_ready$PatientID<-substr(maf_append_purity_ready$Tumor_Sample_Barcode,1,12)
any(colnames(clinical) %in% colnames(maf_append_purity_ready))
# Select ids in clinical file present in maf file
clinical_filter<-subset(clinical, bcr_patient_barcode %in% maf_append_purity_ready$PatientID)
dim(clinical)
dim(clinical_filter)
# Merge files
maf_append_purity_ready_clinical<-merge(maf_append_purity_ready,clinical_filter,by.x="PatientID",by.y="bcr_patient_barcode", all = TRUE)
dim(maf_append_purity_ready_clinical)
dim(maf_append_purity_ready)
head(maf_append_purity_ready_clinical)
# Transform blank spaces to NA!!!
maf_append_purity_ready_clinical[maf_append_purity_ready_clinical==""] <- NA
head(maf_append_purity_ready_clinical)
# Save file
saveRDS(maf_append_purity_ready_clinical,paste(dataDir,"/Resources/maf_append_purity_ready_clinical.rds",sep = "/"))
table(maf_append_purity_ready_clinical$MC3_Overlap) #92.4% (18413928) variants overlap MC3 (MC3_Overlap == TRUE)

################################################################################
### Step 5: Plot Figure 1 
################################################################################

## Read data
maf_append_purity_ready_clinical<-readRDS(paste(dataDir,"/Resources/maf_append_purity_ready_clinical.rds",sep = "/"))
# Put caller column as last column
head(maf_append_purity_ready_clinical)
colnames(maf_append_purity_ready_clinical)
maf_append_purity_ready_clinical<-maf_append_purity_ready_clinical[,c(1:134,137:169,135:136)]
# Prepare for Upset plot
maf_append_purity_ready_clinical$ajcc_stage_parsed<-ifelse(is.na(maf_append_purity_ready_clinical$ajcc_stage_parsed),"Not Available",maf_append_purity_ready_clinical$ajcc_stage_parsed)
table(maf_append_purity_ready_clinical$ajcc_stage_parsed)
sum(is.na(maf_append_purity_ready_clinical$GDC_Validation_Status))

## Check data
dim(maf_append_purity_ready_clinical)
length(which(maf_append_purity_ready_clinical$Reference_Allele == maf_append_purity_ready_clinical$Tumor_Seq_Allele1)) #99.25%
length(which(maf_append_purity_ready_clinical$Tumor_Seq_Allele1 == maf_append_purity_ready_clinical$Tumor_Seq_Allele2)) #0.25%
length(which(maf_append_purity_ready_clinical$Reference_Allele == maf_append_purity_ready_clinical$Tumor_Seq_Allele2)) #0

## Reshape data for plotting
maf_append_purity_ready_clinical_dcast<-reshape2::dcast(maf_append_purity_ready_clinical, name+MC3_Overlap+t_vaf_purity_ploidy_final_mean+total_coverage_mean+ajcc_stage_parsed+GDC_Validation_Status~Caller, fun.aggregate = length)
head(maf_append_purity_ready_clinical_dcast)
colnames(maf_append_purity_ready_clinical_dcast)[c(1:13)]<-c("Variant","MC3","VAF","Coverage","Stage","Validation","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")

table(maf_append_purity_ready_clinical_dcast$MC3)
table(maf_append_purity_ready_clinical_dcast$Stage)
table(maf_append_purity_ready_clinical_dcast$Validation)

callers_maf<-c("MC3","Validation",callers)
maf_append_purity_ready_clinical_dcast$MC3<-ifelse(maf_append_purity_ready_clinical_dcast$MC3 == TRUE, 1, 0)
maf_append_purity_ready_clinical_dcast$Validation<-ifelse(maf_append_purity_ready_clinical_dcast$Validation == "Valid", 1, 0)

maf_append_purity_ready_clinical_dcast[c(callers_maf)]<-maf_append_purity_ready_clinical_dcast[c(callers_maf)] != 0
maf_append_purity_ready_clinical_dcast$MC3<-factor(maf_append_purity_ready_clinical_dcast$MC3,levels = c("TRUE","FALSE"))
maf_append_purity_ready_clinical_dcast$Validation<-factor(maf_append_purity_ready_clinical_dcast$Validation,levels = c("TRUE","FALSE"))
head(maf_append_purity_ready_clinical_dcast)

## Upset plot
# Select coverage thresholds for plotting
summary(maf_append_purity_ready_clinical_dcast$Coverage) #Median Coverage 68 Mean Coverage 103
maf_append_purity_ready_clinical_dcast$Coverage_level<-ifelse(maf_append_purity_ready_clinical_dcast$Coverage<40,"<40",
                                                              ifelse(maf_append_purity_ready_clinical_dcast$Coverage<=70,"40-70",
                                                                     ifelse(maf_append_purity_ready_clinical_dcast$Coverage<=125,">70-125",">125")))
head(maf_append_purity_ready_clinical_dcast)
maf_append_purity_ready_clinical_dcast$Coverage_level<-factor(maf_append_purity_ready_clinical_dcast$Coverage_level,levels = c(">125",">70-125","40-70","<40"))
table(maf_append_purity_ready_clinical_dcast$Coverage_level)
# Prepare stage information por plotting
maf_append_purity_ready_clinical_dcast$Stage<-factor(maf_append_purity_ready_clinical_dcast$Stage,levels = c("Stage I","Stage II","Stage III","Stage IV", "Not Available"))
table(maf_append_purity_ready_clinical_dcast$Stage)
# Save data
saveRDS(maf_append_purity_ready_clinical_dcast,paste(dataDir,"Resources/maf_append_purity_ready_clinical_dcast.rds",sep = "/"))

## Read data
maf_append_purity_ready_clinical_dcast<-readRDS(paste(dataDir,"Resources/maf_append_purity_ready_clinical_dcast.rds",sep = "/"))
head(maf_append_purity_ready_clinical_dcast)
colnames(maf_append_purity_ready_clinical_dcast)
# Check data prior to plotting
maf_append_purity_ready_clinical_dcast$sum<-apply(maf_append_purity_ready_clinical_dcast[,c(9:11,13)],1,function(x) sum(x))
length(which(maf_append_purity_ready_clinical_dcast$Consensus2 == TRUE & maf_append_purity_ready_clinical_dcast$sum <2))
length(which(maf_append_purity_ready_clinical_dcast$Consensus2 == FALSE & maf_append_purity_ready_clinical_dcast$sum >1))
length(which(maf_append_purity_ready_clinical_dcast$Consensus3 == TRUE & maf_append_purity_ready_clinical_dcast$sum <3))
length(which(maf_append_purity_ready_clinical_dcast$Consensus3 == FALSE & maf_append_purity_ready_clinical_dcast$sum >2))
length(which(maf_append_purity_ready_clinical_dcast$Union == TRUE & maf_append_purity_ready_clinical_dcast$sum <1))

## PLOT FIGURE 1
# Set color palette
cbPalette <- c("#E69F00","#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#  Upset plot
Figure_maf<-
  (
    upset(maf_append_purity_ready_clinical_dcast, callers, name=NULL, width_ratio=0.2, height_ratio = 0.8,  n_intersections = 14, keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("VAF"=(ggplot(mapping=aes(y=VAF))+geom_jitter(aes(color=Coverage_level),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(0.85,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>125'="cadetblue1", '>70-125'="aquamarine1", '40-70'="antiquewhite3", '<40'="antiquewhite")))),
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),text=list(size=2.25),mapping=aes(fill=MC3)) + scale_fill_manual(name="MC3 Overlap",values=c('TRUE'="darkgoldenrod1", 'FALSE'="darkblue")) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=MC3)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 1.85)+ expand_limits(y=4500000)+theme(legend.position = "none")+scale_fill_manual(values=c('TRUE'="darkgoldenrod1", 'FALSE'="darkblue"))))
    + ggtitle("Variant Calling Results") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(plot = Figure_maf, width = 9, height = 5, dpi = 300, filename = paste(dataDir,"Figures/Figure1.png",sep = "/"))

# Save data
saveRDS(Figure_maf,paste(dataDir,"Figures/Figure1.rds",sep="/"))
#===============================================================================