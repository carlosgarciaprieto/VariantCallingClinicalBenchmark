#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate data for Figure 5, Figure S6 and Figure S7
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
### Step 1: Read MOAlmanac results
################################################################################

## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          "LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          "THCA","THYM","UCEC","UCS","UVM")

caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")

# We select all files with actionable variants
path<-paste(dataDir,"MOAlmanac",sep = "/")
files<-list.files(path=path, pattern = "*actionable.txt", full.names = TRUE)
length(files)
# Read data
moalmanac_final<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    moalmanac<-NULL
    moalmanac<-as.data.frame(fread(paste0(dataDir,"/Moalmanac/",cancer[i],"_",caller[k],".actionable.txt")))
    moalmanac$Cancer<-cancer[i]
    moalmanac$Caller<-caller[k]
    moalmanac_final<-rbind(moalmanac_final,moalmanac)
  }
}
head(moalmanac_final)
# Save data
saveRDS(moalmanac_final,paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable.rds",sep = "/"))

################################################################################
### Step 2: Prepare MOAlmanac data
################################################################################

# Read data
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable.rds",sep = "/"))
# Compute tumor_f (variant allelic fraction) and total_coverage means by feature display
tumor_f_mean<-aggregate(moalmanac_final$tumor_f, list(moalmanac_final$feature_display), FUN=mean)
total_coverage_mean<-aggregate(moalmanac_final$total_coverage, list(moalmanac_final$feature_display), FUN=mean)
colnames(tumor_f_mean)<-c("feature_display","tumor_f_mean")
colnames(total_coverage_mean)<-c("feature_display","total_coverage_mean")
total_coverage_mean$total_coverage_mean<-round(total_coverage_mean$total_coverage_mean)
head(tumor_f_mean)
head(total_coverage_mean)
# Merge both means
tumor_f_coverage_mean<-merge(tumor_f_mean,total_coverage_mean,by="feature_display")
head(tumor_f_coverage_mean)
dim(tumor_f_coverage_mean)
# Check dimensions
dim(tumor_f_mean)
dim(total_coverage_mean)
length(unique(moalmanac_final$feature_display))
# Merge with MOAlmanac results
moalmanac_final_ready<-merge(moalmanac_final,tumor_f_coverage_mean, by = "feature_display")
dim(moalmanac_final_ready)
dim(moalmanac_final)
# Transform blank spaces to NA!!!
moalmanac_final_ready[moalmanac_final_ready==""] <- NA
# Save data
saveRDS(moalmanac_final_ready,paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready.rds",sep = "/"))

################################################################################
### Step 3: Manually curate MOAlmanac results to annotate cancer type for each feature display
################################################################################

## Therapeutic sensitivity
# Read data
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready.rds",sep = "/"))
# Check all sensitive descriptions
sensitive_descriptions<-unique(moalmanac_final$sensitive_description)
# Manually curate and annotate cancer type
moalmanac_final$sensitive_description_cancer<-case_when(moalmanac_final$sensitive_description==sensitive_descriptions[2] ~ "PANCANCER"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[3] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[4] ~ "LUAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[5] ~ "OV;Other"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[6] ~ "Solid tumor"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[7] ~ "PRAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[8] ~ "COAD;READ"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[9] ~ "SKCM"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[10] ~ "PAAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[11] ~ "OV;Ovarian Epithelial Tumor;High-Grade Serous Fallopian Tube Cancer;Peritoneal Serous Carcinoma"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[12] ~ "PAAD;Pancreatic cancers"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[13] ~ "OV;Ovarian Epithelial Tumor;High-Grade Serous Fallopian Tube Cancer;Peritoneal Serous Carcinoma"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[14] ~ "SKCM"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[15] ~ "PANCANCER;COAD;READ"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[16] ~ "BLCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[17] ~ "STAD;Gastric"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[18] ~ "PAAD;Pancreatic cancers"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[19] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[20] ~ "BLCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[21] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[22] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[23] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[24] ~ "BRCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[25] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[26] ~ "BLCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[27] ~ "CHOL"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[28] ~ "Follicular lymphoma"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[29] ~ "BRCA;Breast Invasive Ductal Carcinoma"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[30] ~ "BLCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[31] ~ "BLCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[32] ~ "LAML"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[33] ~ "BRCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[34] ~ "LAML"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[35] ~ "LAML"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[36] ~ "ALL;Acute Lymphoblastic Leukemia"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[37] ~ "ALL;Acute Lymphoblastic Leukemia"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[38] ~ "GIST;Gastrointestinal Stromal Tumor"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[39] ~ "HNSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[40] ~ "GIST;Gastrointestinal Stromal Tumor"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[41] ~ "SKCM"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[42] ~ "OV"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[43] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[44] ~ "LUAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[45] ~ "HNSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[46] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[47] ~ "HNSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[48] ~ "PANCANCER;COAD;READ"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[49] ~ "Myeloproliferative neoplasm"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[50] ~ "SKCM;PANCANCER"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[51] ~ "KIRC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[52] ~ "GIST;Gastrointestinal Stromal Tumor"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[53] ~ "BRCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[54] ~ "BRCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[55] ~ "PANCANCER;PRAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[56] ~ "PANCANCER;COAD;READ"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[57] ~ "PANCANCER;PRAD"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[58] ~ "Medullary Thyroid Cancer;THCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[59] ~ "Medullary Thyroid Cancer;THCA"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[60] ~ "Non-small cell lung cancer;LUAD;LUSC"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[61] ~ "Myelodysplastic Syndromes;MDS"
                                                        ,moalmanac_final$sensitive_description==sensitive_descriptions[62] ~ "Subependymal Giant Cell Astrocytoma;SEGA;LGG"
)


View(table(moalmanac_final$sensitive_description,moalmanac_final$sensitive_description_cancer))
anyNA(moalmanac_final$sensitive_description_cancer)

## Therapeutic resistance
# Check all resistance descriptions
resistance_descriptions<-unique(moalmanac_final$resistance_description)
# Manually curate and annotate cancer type
moalmanac_final$resistance_description_cancer<-ifelse(moalmanac_final$resistance_description==resistance_descriptions[1],"Chronic Myelogenous Leukemia;CML;LCML",
                                                      ifelse(moalmanac_final$resistance_description==resistance_descriptions[3],"Non-small cell lung cancer;LUAD;LUSC",
                                                             ifelse(moalmanac_final$resistance_description==resistance_descriptions[4],"PRAD",
                                                                    ifelse(moalmanac_final$resistance_description==resistance_descriptions[5],"Chronic Myelogenous Leukemia;CML;LCML",
                                                                           ifelse(moalmanac_final$resistance_description==resistance_descriptions[6],"COAD;READ",
                                                                                  ifelse(moalmanac_final$resistance_description==resistance_descriptions[7],"OV",
                                                                                         ifelse(moalmanac_final$resistance_description==resistance_descriptions[8],"OV",
                                                                                                ifelse(moalmanac_final$resistance_description==resistance_descriptions[9],"Non-small cell lung cancer;LUAD;LUSC",
                                                                                                       ifelse(moalmanac_final$resistance_description==resistance_descriptions[10],"Non-small cell lung cancer;LUAD;LUSC",
                                                                                                              ifelse(moalmanac_final$resistance_description==resistance_descriptions[11],"BRCA",
                                                                                                                     ifelse(moalmanac_final$resistance_description==resistance_descriptions[12],"LAML;AML",
                                                                                                                            ifelse(moalmanac_final$resistance_description==resistance_descriptions[13],"READ",
                                                                                                                                   ifelse(moalmanac_final$resistance_description==resistance_descriptions[14],"Solid tumor",
                                                                                                                                          ifelse(moalmanac_final$resistance_description==resistance_descriptions[15],"COAD;READ",
                                                                                                                                                 ifelse(moalmanac_final$resistance_description==resistance_descriptions[16],"SKCM",
                                                                                                                                                        ifelse(moalmanac_final$resistance_description==resistance_descriptions[17],"SKCM",
                                                                                                                                                               ifelse(moalmanac_final$resistance_description==resistance_descriptions[18],"SKCM",
                                                                                                                                                                      ifelse(moalmanac_final$resistance_description==resistance_descriptions[19],"SKCM",
                                                                                                                                                                             ifelse(moalmanac_final$resistance_description==resistance_descriptions[20],"LUSC",
                                                                                                                                                                                    ifelse(moalmanac_final$resistance_description==resistance_descriptions[21],"Uterine Leiomyoma;ULM;",
                                                                                                                                                                                           ifelse(moalmanac_final$resistance_description==resistance_descriptions[22],"Medullary Thyroid Cancer;THCA",
                                                                                                                                                                                                  ifelse(moalmanac_final$resistance_description==resistance_descriptions[23],"Myelodysplastic Syndromes;MDS"
                                                                                                                                                                                                         ,NA))))))))))))))))))))))


## Therapeutic prognosis
# Check all prognostic descriptions
prognostic_descriptions<-unique(moalmanac_final$prognostic_description)
# Manually curate and annotate cancer type
moalmanac_final$prognostic_description_cancer<-ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[2],"PAAD",
                                                      ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[3],"Chronic Myelomonocytic Leukemia;CMML",
                                                             ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[4],"OV;Ovarian Epithelial Tumor;Fallopian Tube Cancer;Peritoneal Serous Carcinoma",
                                                                    ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[5],"Chronic Myelomonocytic Leukemia;CMML",
                                                                           ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[6],"COAD;READ",
                                                                                  ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[7],"COAD;READ",
                                                                                         ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[8],"Medulloblastoma;MBL",
                                                                                                ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[9],"Myelodysplastic Syndromes;MDS",
                                                                                                       ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[10],"May be familial in rare cases",
                                                                                                              ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[11],"Glioma;LGG;GBM",
                                                                                                                     ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[12],"LAML;AML",
                                                                                                                            ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[13],"Glioma;LGG;GBM",
                                                                                                                                   ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[14],"ALL;Acute Lymphoblastic Leukemia",
                                                                                                                                          ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[15],"Myeloproliferative neoplasm",
                                                                                                                                                 ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[16],"HNSC",
                                                                                                                                                        ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[17],"Myelodysplastic Syndromes;MDS",
                                                                                                                                                               ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[18],"PAAD",
                                                                                                                                                                      ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[19],"SKCM",
                                                                                                                                                                             ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[20],"Chronic Myelomonocytic Leukemia;CMML",
                                                                                                                                                                                    ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[21],"PRAD",
                                                                                                                                                                                           ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[22],"COAD;READ",
                                                                                                                                                                                                  ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[23],"Myelodysplastic Syndromes;MDS",
                                                                                                                                                                                                         ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[24],"PAAD",
                                                                                                                                                                                                                ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[25],"Chronic Myelomonocytic Leukemia;CML",
                                                                                                                                                                                                                       ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[26],"Myelodysplastic Syndromes;MDS",
                                                                                                                                                                                                                              ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[27],"Myelodysplastic Syndromes;MDS",
                                                                                                                                                                                                                                     ifelse(moalmanac_final$prognostic_description==prognostic_descriptions[28],"Myelodysplastic Syndromes;MDS",
                                                                                                                                                                                                                                            NA)))))))))))))))))))))))))))


#Save results
saveRDS(moalmanac_final,paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready_curated.rds",sep = "/"))

################################################################################
### Step 4: Add purity/ploidy data to MOAlmanac results
################################################################################

# Read MOAlmanac ready to use results
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready_curated.rds",sep = "/"))
head(moalmanac_final)
# Read adjusted VAF data
maf_append_purity_ready_clinical<-readRDS(paste(dataDir,"/Resources/maf_append_purity_ready_clinical.rds",sep = "/"))
head(maf_append_purity_ready_clinical)
## Merge MOAlmanac results with adjusted VAF data
# Select relevant purity/ploidy information from samples present in MOAlmanac results file
maf_append_purity_ready_clinical_filter<-maf_append_purity_ready_clinical[maf_append_purity_ready_clinical$Tumor_Sample_Barcode %in% moalmanac_final$tumor_sample_barcode,c("Tumor_Sample_Barcode","ajcc_pathologic_tumor_stage","ajcc_stage_parsed","purity","ploidy","Cancer DNA fraction","Subclonal genome fraction","Genome doublings")]
maf_append_purity_ready_clinical_filter<-maf_append_purity_ready_clinical_filter[!duplicated(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode)]
dim(maf_append_purity_ready_clinical_filter)
head(maf_append_purity_ready_clinical_filter)
length(unique(moalmanac_final$tumor_sample_barcode))
sum(is.na(moalmanac_final$tumor_sample_barcode))
length(unique(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode))
sum(is.na(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode))
# Merge
any(colnames(maf_append_purity_ready_clinical_filter) %in% colnames(moalmanac_final))
length(which(unique(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode) %in% moalmanac_final$tumor_sample_barcode))
moalmanac_final_adjust<-merge(moalmanac_final,maf_append_purity_ready_clinical_filter, by.x = "tumor_sample_barcode",by.y = "Tumor_Sample_Barcode", all = TRUE)
dim(moalmanac_final_adjust)
dim(moalmanac_final)
head(moalmanac_final_adjust)
#Save results
saveRDS(moalmanac_final_adjust,paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready_curated_adjust.rds",sep = "/"))

################################################################################
### Step 5: Add annotation (MC3 overlap and validation status) to MOAlmanac results
################################################################################

# Read MOAlmanac with purity/ploidy data
moalmanac_final<-readRDS(paste(dataDir,"moalmanac_actionable_ready_curated_adjust.rds",sep = "/"))
head(moalmanac_final)
head(maf_append_purity_ready_clinical_filter)
## Merge MOAlmanac results with annotation data
# Select  samples in MOAlmanac results file
maf_append_purity_ready_clinical_annot<-maf_append_purity_ready_clinical[maf_append_purity_ready_clinical$Tumor_Sample_Barcode %in% moalmanac_final$tumor_sample_barcode,c("Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","MC3_Overlap","GDC_Validation_Status")]
# Rename to same nomenclature
table(moalmanac_final$alteration_type)
levels(maf_append_purity_ready_clinical_annot$Variant_Classification)
table(maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("In_Frame_Del","Deletion",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Frame_Shift_Ins","Frameshift",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Frame_Shift_Del","Frameshift",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Missense_Mutation","Missense",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Nonsense_Mutation","Nonsense",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Splice_Site","Splice Site",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("In_Frame_Ins","Insertion",maf_append_purity_ready_clinical_annot$Variant_Classification)
maf_append_purity_ready_clinical_annot$Variant_Classification<-gsub("Nonstop_Mutation","Nonstop",maf_append_purity_ready_clinical_annot$Variant_Classification)
table(maf_append_purity_ready_clinical_annot$Variant_Classification)
# Create key to merge (merge by feature and sample)
head(maf_append_purity_ready_clinical_annot)
head(moalmanac_final)
moalmanac_final$MergeID<-paste(moalmanac_final$feature_display,moalmanac_final$tumor_sample_barcode,sep = ":")
maf_append_purity_ready_clinical_annot$MergeID<-paste0(maf_append_purity_ready_clinical_annot$Hugo_Symbol," ",maf_append_purity_ready_clinical_annot$HGVSp_Short," ","(",maf_append_purity_ready_clinical_annot$Variant_Classification,")",":",maf_append_purity_ready_clinical_annot$Tumor_Sample_Barcode)
# Remove duplicates
maf_append_purity_ready_clinical_annot<-maf_append_purity_ready_clinical_annot[!duplicated(maf_append_purity_ready_clinical_annot$MergeID)]
dim(maf_append_purity_ready_clinical_annot)
head(maf_append_purity_ready_clinical_annot)
# Select only variants in MOAlmanac results
maf_append_purity_ready_clinical_annot_filter<-maf_append_purity_ready_clinical_annot[maf_append_purity_ready_clinical_annot$MergeID %in% moalmanac_final$MergeID,]
dim(maf_append_purity_ready_clinical_annot_filter)
length(unique(moalmanac_final$MergeID))
# Finally merge both files
moalmanac_final_annot<-merge(moalmanac_final,maf_append_purity_ready_clinical_annot_filter,by="MergeID",all = TRUE)
dim(moalmanac_final_annot)
dim(moalmanac_final)
# Remove blank spaces
moalmanac_final_annot[moalmanac_final_annot==""] <- NA
# Adjust purity
sum(is.na(moalmanac_final_annot$`Cancer DNA fraction`))/(dim(moalmanac_final_annot)[1])*100 #35232 variants 16.5% variants missing Cancer DNA fraction
moalmanac_final_annot$t_vaf_purity<-moalmanac_final_annot$tumor_f/moalmanac_final_annot$`Cancer DNA fraction`
length(which(moalmanac_final_annot$t_vaf_purity>1)) #9716 variants
length(which(moalmanac_final_annot$t_vaf_purity>1))/(dim(moalmanac_final_annot)[1])*100 #4.5% of variants
length(which(moalmanac_final_annot$t_vaf_purity<1)) #168408 variants
summary(moalmanac_final_annot$t_vaf_purity)
summary(moalmanac_final_annot$tumor_f)
dim(moalmanac_final_annot)
# Adjust ploidy
moalmanac_final_annot$t_vaf_purity_ploidy<-moalmanac_final_annot$t_vaf_purity/(moalmanac_final_annot$ploidy/2)
length(which(round(moalmanac_final_annot$t_vaf_purity_ploidy,1)>1)) #4167 variants
length(which(round(moalmanac_final_annot$t_vaf_purity_ploidy,1)>1))/(dim(moalmanac_final_annot)[1])*100 #1.95% of variants
length(which(moalmanac_final_annot$t_vaf_purity_ploidy<1)) #172187 variants
summary(moalmanac_final_annot$t_vaf_purity_ploidy)
summary(moalmanac_final_annot$t_vaf_purity)
summary(moalmanac_final_annot$t_vaf)
dim(moalmanac_final_annot)
# Finally replace adjusted VAFs > 1 (samples are matched by plate, no by exact matches) and NAs with unadjusted VAFs
sum(is.na(moalmanac_final_annot$t_vaf_purity_ploidy))/(dim(moalmanac_final_annot)[1])*100 #16.5% NAs
moalmanac_final_annot$t_vaf_purity_ploidy_final<-ifelse(is.na(moalmanac_final_annot$t_vaf_purity_ploidy), moalmanac_final_annot$tumor_f,ifelse(round(moalmanac_final_annot$t_vaf_purity_ploidy,1) >1, moalmanac_final_annot$tumor_f,moalmanac_final_annot$t_vaf_purity_ploidy))
dim(moalmanac_final_annot)
sum(is.na(moalmanac_final_annot$tumor_f)) #560 samples missing tumor_f info
head(moalmanac_final_annot)
##Save file with purity ploidy data
saveRDS(moalmanac_final_annot,paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready_curated_adjust_annot.rds",sep = "/"))

################################################################################
### Step 6: Reshape data (Aggregate VAF and coverage means by variant)
################################################################################

# Read MOAlmanac with purity/ploidy and annotation data
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_actionable_ready_curated_adjust_annot.rds",sep = "/"))
head(moalmanac_final)
dim(moalmanac_final)
#Compute VAF and coverage means by variant
moalmanac_tumor_f_mean<-aggregate(moalmanac_final$t_vaf_purity_ploidy_final, list(moalmanac_final$MergeID), FUN=mean)
moalmanac_total_coverage_mean<-aggregate(moalmanac_final$total_coverage, list(moalmanac_final$MergeID), FUN=mean)
colnames(moalmanac_tumor_f_mean)<-c("MergeID","tumor_f_mean")
colnames(moalmanac_total_coverage_mean)<-c("MergeID","total_coverage_mean")
moalmanac_total_coverage_mean$total_coverage_mean<-round(moalmanac_total_coverage_mean$total_coverage_mean)
head(moalmanac_tumor_f_mean)
dim(moalmanac_tumor_f_mean)
length(unique(moalmanac_final$MergeID))
head(moalmanac_total_coverage_mean)
dim(moalmanac_total_coverage_mean)
#Merge both tumor_f (variant allelic fraction) and coverage means data
moalmanac_tumor_f_coverage_mean<-merge(moalmanac_tumor_f_mean,moalmanac_total_coverage_mean,by="MergeID")
head(moalmanac_tumor_f_coverage_mean)
dim(moalmanac_tumor_f_coverage_mean)
# Save results
saveRDS(moalmanac_tumor_f_coverage_mean,paste(dataDir,"Resources/MOAlmanac/moalmanac_tumor_f_coverage_mean.rds",sep = "/"))

## Merge with MOAlmanac file
# Remove previous unadjusted means
moalmanac<-merge(moalmanac_final[,-c(45,46)],moalmanac_tumor_f_coverage_mean, by = "MergeID")
dim(moalmanac)
dim(moalmanac_final)
#Transform blank spaces to NA!!!
moalmanac[moalmanac==""] <- NA
head(moalmanac)
#Save file
saveRDS(moalmanac,paste(dataDir,"/Resources/MOAlmanac/moalmanac_results_adjusted_annotated.rds",sep = "/"))
#===============================================================================