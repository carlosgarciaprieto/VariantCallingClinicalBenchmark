#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure S8
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
### Step 1: Read MOAlmanac MSI results
################################################################################

## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          "ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          "LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          "THCA","THYM","UCEC","UCS","UVM")

caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")

# We select all files with MSI actionable variants
path<-paste(dataDir,"MOAlmanac",sep = "/")
files<-list.files(path=path, pattern = "*msi_variants.txt", full.names = TRUE)
length(files)
# Read data
moalmanac_msi_final<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    moalmanac_msi<-NULL
    moalmanac_msi<-as.data.frame(fread(paste0(dataDir,"/Moalmanac/",cancer[i],"_",caller[k],".msi_variants.txt")))
    if (nrow(moalmanac_msi)>0){
      moalmanac_msi$Cancer<-cancer[i]
      moalmanac_msi$Caller<-caller[k]
      moalmanac_msi_final<-rbind(moalmanac_msi_final,moalmanac_msi)
    }
    else {print(paste(i,cancer[i],caller[k],"No MSI variants detected",sep = "_"))}
  }
}
head(moalmanac_msi_final)
# Add key ID column
moalmanac_msi_final$feature_display<-paste0(moalmanac_msi_final$feature," ",moalmanac_msi_final$alteration," ","(",moalmanac_msi_final$alteration_type,")")
# Save data
saveRDS(moalmanac_msi_final,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi.rds",sep = "/"))

################################################################################
### Step 2: Prepare MOAlmanac MSI data
################################################################################

# Read MSI results
moalmanac_msi_final<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi.rds",sep = "/"))
head(moalmanac_msi_final)
# Read adjusted VAF data
maf_append_purity_ready_clinical<-readRDS(paste(dataDir,"/Resources/maf_append_purity_ready_clinical.rds",sep = "/"))
head(maf_append_purity_ready_clinical)

## Merge MOAlmanac MIS results with adjusted VAF data
# Select samples in MOAlmanac MSI results file
maf_append_purity_ready_clinical_filter<-maf_append_purity_ready_clinical[maf_append_purity_ready_clinical$Tumor_Sample_Barcode %in% moalmanac_msi_final$tumor_sample_barcode,c("Tumor_Sample_Barcode","ajcc_pathologic_tumor_stage","ajcc_stage_parsed","purity","ploidy","Cancer DNA fraction","Subclonal genome fraction","Genome doublings")]
maf_append_purity_ready_clinical_filter<-maf_append_purity_ready_clinical_filter[!duplicated(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode)]
dim(maf_append_purity_ready_clinical_filter)
head(maf_append_purity_ready_clinical_filter)
length(unique(moalmanac_msi_final$tumor_sample_barcode))
sum(is.na(moalmanac_msi_final$tumor_sample_barcode))
length(unique(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode))
sum(is.na(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode))
# Merge both files
any(colnames(maf_append_purity_ready_clinical_filter) %in% colnames(moalmanac_msi_final))
length(which(unique(maf_append_purity_ready_clinical_filter$Tumor_Sample_Barcode) %in% moalmanac_msi_final$tumor_sample_barcode))
moalmanac_msi_final_adjust<-merge(moalmanac_msi_final,maf_append_purity_ready_clinical_filter, by.x = "tumor_sample_barcode",by.y = "Tumor_Sample_Barcode", all = TRUE)
dim(moalmanac_msi_final_adjust)
dim(moalmanac_msi_final)
head(moalmanac_msi_final_adjust)
# Save results
saveRDS(moalmanac_msi_final_adjust,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_ready_curated_adjust.rds",sep = "/"))

################################################################################
### Step 3: Add annotation (MC3 overlap and validation status) and purity/ploidy data to MOAlmanac MSI results
################################################################################

## Add annotation to MOAlmanac MSI results (MC3 overlap and validation status) 
# Read MOAlmanac with purity/ploidy data
moalmanac_msi_final<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_ready_curated_adjust.rds",sep = "/"))
head(moalmanac_msi_final)
head(maf_append_purity_ready_clinical_filter)

## Merge MOAlmanac MSI results with annotation data
# Select samples in MOAlmanac MSI results file
maf_append_purity_ready_clinical_annot<-maf_append_purity_ready_clinical[maf_append_purity_ready_clinical$Tumor_Sample_Barcode %in% moalmanac_msi_final$tumor_sample_barcode,c("Tumor_Sample_Barcode","Hugo_Symbol","HGVSp_Short","Variant_Classification","MC3_Overlap","GDC_Validation_Status")]
# Rename to same nomenclature
table(moalmanac_msi_final$alteration_type)
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
# Create key ID column to merge
head(maf_append_purity_ready_clinical_annot)
head(moalmanac_msi_final)
moalmanac_msi_final$MergeID<-paste(moalmanac_msi_final$feature_display,moalmanac_msi_final$tumor_sample_barcode,sep = ":")
maf_append_purity_ready_clinical_annot$MergeID<-paste0(maf_append_purity_ready_clinical_annot$Hugo_Symbol," ",maf_append_purity_ready_clinical_annot$HGVSp_Short," ","(",maf_append_purity_ready_clinical_annot$Variant_Classification,")",":",maf_append_purity_ready_clinical_annot$Tumor_Sample_Barcode)
# Remove duplicates
maf_append_purity_ready_clinical_annot<-maf_append_purity_ready_clinical_annot[!duplicated(maf_append_purity_ready_clinical_annot$MergeID)]
dim(maf_append_purity_ready_clinical_annot)
head(maf_append_purity_ready_clinical_annot)
# Select only variants in MOAlmanac MSI results
maf_append_purity_ready_clinical_annot_filter<-maf_append_purity_ready_clinical_annot[maf_append_purity_ready_clinical_annot$MergeID %in% moalmanac_msi_final$MergeID,]
dim(maf_append_purity_ready_clinical_annot_filter)
length(unique(moalmanac_msi_final$MergeID))

## Finally merge both files
moalmanac_msi_final_annot<-merge(moalmanac_msi_final,maf_append_purity_ready_clinical_annot_filter,by="MergeID",all = TRUE)
dim(moalmanac_msi_final_annot)
dim(moalmanac_msi_final)
# Remove blank spaces
moalmanac_msi_final_annot[moalmanac_msi_final_annot==""] <- NA
# Adjust purity
sum(is.na(moalmanac_msi_final_annot$`Cancer DNA fraction`))/(dim(moalmanac_msi_final_annot)[1])*100 #792 variants 14.8% variants missing Cancer DNA fraction
moalmanac_msi_final_annot$t_vaf_purity<-moalmanac_msi_final_annot$tumor_f/moalmanac_msi_final_annot$`Cancer DNA fraction`
length(which(moalmanac_msi_final_annot$t_vaf_purity>1)) #222 variants
length(which(moalmanac_msi_final_annot$t_vaf_purity>1))/(dim(moalmanac_msi_final_annot)[1])*100 #4.15% of variants
length(which(moalmanac_msi_final_annot$t_vaf_purity<1)) #4327 variants
summary(moalmanac_msi_final_annot$t_vaf_purity)
summary(moalmanac_msi_final_annot$tumor_f)
dim(moalmanac_msi_final_annot)
#Adjust ploidy
moalmanac_msi_final_annot$t_vaf_purity_ploidy<-moalmanac_msi_final_annot$t_vaf_purity/(moalmanac_msi_final_annot$ploidy/2)
length(which(round(moalmanac_msi_final_annot$t_vaf_purity_ploidy,1)>1)) #142 variants
length(which(round(moalmanac_msi_final_annot$t_vaf_purity_ploidy,1)>1))/(dim(moalmanac_msi_final_annot)[1])*100 #2.66% of variants
length(which(moalmanac_msi_final_annot$t_vaf_purity_ploidy<1)) #4325 variants
summary(moalmanac_msi_final_annot$t_vaf_purity_ploidy)
summary(moalmanac_msi_final_annot$t_vaf_purity)
summary(moalmanac_msi_final_annot$tumor_f)
dim(moalmanac_msi_final_annot)
#Finally replace VAFs > 1 (samples are matched by plate, no exact matches) and NAs with unadjusted VAFs
sum(is.na(moalmanac_msi_final_annot$t_vaf_purity_ploidy))/(dim(moalmanac_msi_final_annot)[1])*100 #14.8% NAs
moalmanac_msi_final_annot$t_vaf_purity_ploidy_final<-ifelse(is.na(moalmanac_msi_final_annot$t_vaf_purity_ploidy), moalmanac_msi_final_annot$tumor_f,ifelse(round(moalmanac_msi_final_annot$t_vaf_purity_ploidy,1) >1, moalmanac_msi_final_annot$tumor_f,moalmanac_msi_final_annot$t_vaf_purity_ploidy))
dim(moalmanac_msi_final_annot)
sum(is.na(moalmanac_msi_final_annot$tumor_f)) #0 samples missing tumor_f (variant allelic fraction) information
head(moalmanac_msi_final_annot)
##Save file with purity ploidy data
saveRDS(moalmanac_msi_final_annot,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_ready_curated_adjust_annot.rds",sep = "/"))

################################################################################
### Step 4: Reshape data (Aggregate VAF and coverage means by variant)
################################################################################

#Read MOAlmanac MSI with purity/ploidy and annotation data
moalmanac_msi_final<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_ready_curated_adjust_annot.rds",sep = "/"))
head(moalmanac_msi_final)
dim(moalmanac_msi_final)
# Compute VAF and coverage means by variant
moalmanac_tumor_f_mean<-aggregate(moalmanac_msi_final$t_vaf_purity_ploidy_final, list(moalmanac_msi_final$MergeID), FUN=mean)
moalmanac_total_coverage_mean<-aggregate(moalmanac_msi_final$total_coverage, list(moalmanac_msi_final$MergeID), FUN=mean)
colnames(moalmanac_tumor_f_mean)<-c("MergeID","tumor_f_mean")
colnames(moalmanac_total_coverage_mean)<-c("MergeID","total_coverage_mean")
moalmanac_total_coverage_mean$total_coverage_mean<-round(moalmanac_total_coverage_mean$total_coverage_mean)
head(moalmanac_tumor_f_mean)
dim(moalmanac_tumor_f_mean)
length(unique(moalmanac_msi_final$MergeID))
head(moalmanac_total_coverage_mean)
dim(moalmanac_total_coverage_mean)
# Merge both tumor_f (variant allelic fraction) and coverage means data
moalmanac_tumor_f_coverage_mean<-merge(moalmanac_tumor_f_mean,moalmanac_total_coverage_mean,by="MergeID")
head(moalmanac_tumor_f_coverage_mean)
dim(moalmanac_tumor_f_coverage_mean)
# Save results
saveRDS(moalmanac_tumor_f_coverage_mean,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_tumor_f_coverage_mean.rds",sep = "/"))

## Merge with MOAlmanac MSI file
moalmanac_msi<-merge(moalmanac_msi_final,moalmanac_tumor_f_coverage_mean, by = "MergeID")
dim(moalmanac_msi)
dim(moalmanac_msi_final)
# Transform blank spaces to NA!!!
moalmanac_msi[moalmanac_msi==""] <- NA
head(moalmanac_msi)
# Save file
saveRDS(moalmanac_msi,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_results_adjusted_annotated.rds",sep = "/"))

################################################################################
### Step 5: Plot Figure S8A
################################################################################

# Read MOAlmanac MSI results with purity/ploidy adjusted and aggregated by variant and sample and annotation data included
moalmanac_msi_final<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_results_adjusted_annotated.rds",sep = "/"))
head(moalmanac_msi_final)
# Remove patient id, tumor_f,total coverage,t_vaf_purity,t_vaf_purity_ploidy and t_vaf_purity_ploidy_final columns for Upset plot
colnames(moalmanac_msi_final)
moalmanac_msi_final<-moalmanac_msi_final[,-c(31,17:18,56:58)]
head(moalmanac_msi_final)
## Prepare data
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
Callers<-toupper(callers)
# Put caller column as last column
head(moalmanac_msi_final)
colnames(moalmanac_msi_final)
moalmanac_msi_final<-moalmanac_msi_final[,c(1:36,39:54,37:38)]
colnames(moalmanac_msi_final)[42:44]<-c("Cancer_DNA_fraction","Subclonal_genome_fraction","Genome_doublings")
# Transform blank spaces to NA!!!
moalmanac_msi_final[moalmanac_msi_final==""] <- NA
# Reshape data
moalmanac_msi_dcast<-reshape2::dcast(moalmanac_msi_final, ...~Caller, fun.aggregate = length)
colnames(moalmanac_msi_dcast)[c(3,54:60)]<-c("Bin","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers_moalmanac<-c("MC3_Overlap","GDC_Validation_Status",callers)
moalmanac_msi_dcast$MC3_Overlap<-ifelse(moalmanac_msi_dcast$MC3_Overlap == TRUE, 1, 0) #80.9% msi actionable variants in MC3
moalmanac_msi_dcast$GDC_Validation_Status<-ifelse(moalmanac_msi_dcast$GDC_Validation_Status == "Valid", 1, 0) 
moalmanac_msi_dcast[callers_moalmanac]<-moalmanac_msi_dcast[callers_moalmanac] != 0
head(moalmanac_msi_dcast)
moalmanac_msi_dcast$MC3_Overlap<-factor(moalmanac_msi_dcast$MC3_Overlap,levels = c("TRUE","FALSE"))
moalmanac_msi_dcast$GDC_Validation_Status<-factor(moalmanac_msi_dcast$GDC_Validation_Status,levels = c("TRUE","FALSE"))
# Save file
saveRDS(moalmanac_msi_dcast,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_dcast.rds",sep = "/"))

## Plot Figure S8A
# Read data
moalmanac_msi_dcast<-readRDS(paste(dataDir,"/Resources/MOAlmanac_MSI/moalmanac_msi_dcast.rds",sep = "/"))
# Set color palette
cbPalette <- c("#E69F00","#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Prepare coverage information for plot
summary(moalmanac_msi_dcast$total_coverage_mean) #Median Coverage 73 Mean Coverage 104
moalmanac_msi_dcast$total_coverage_mean_factor<-ifelse(moalmanac_msi_dcast$total_coverage_mean<45,"<45",
                                                       ifelse(moalmanac_msi_dcast$total_coverage_mean<=75,"45-75",
                                                              ifelse(moalmanac_msi_dcast$total_coverage_mean<=130,">75-130",">130")))
head(moalmanac_msi_dcast)
moalmanac_msi_dcast$total_coverage_mean_factor<-factor(moalmanac_msi_dcast$total_coverage_mean_factor,levels = c(">130",">75-130","45-75","<45"))
table(moalmanac_msi_dcast$total_coverage_mean_factor)
# Format stage information
moalmanac_msi_dcast$ajcc_stage_parsed<-factor(moalmanac_msi_dcast$ajcc_stage_parsed,levels = c("Stage I","Stage II","Stage III","Stage IV", "Not Available"))
table(moalmanac_msi_dcast$ajcc_stage_parsed)
# Format Bin column
moalmanac_msi_dcast$Bin_format<-ifelse(moalmanac_msi_dcast$almanac_bin>0 & moalmanac_msi_dcast$Bin=="Biologically Relevant","Biologically Relevant",ifelse(moalmanac_msi_dcast$almanac_bin>0 & moalmanac_msi_dcast$Bin=="Investigate Actionability","Investigate Actionability",ifelse(moalmanac_msi_dcast$cancerhotspots_bin>0,"Cancer Hotspot",ifelse(moalmanac_msi_dcast$cancerhotspots3D_bin>0,"Cancer Hotspot 3D",ifelse(moalmanac_msi_dcast$cgc_bin>0,"Cancer Gene Census",ifelse(moalmanac_msi_dcast$gsea_pathways_bin>0,"Cancer Pathway",ifelse(moalmanac_msi_dcast$gsea_modules_bin>0,"Cancer Module",ifelse(moalmanac_msi_dcast$cosmic_bin>0,"Cosmic",NA))))))))
table(moalmanac_msi_dcast$Bin_format)
sum(is.na(moalmanac_msi_dcast$Bin_format))
moalmanac_msi_dcast$Bin_format<-factor(moalmanac_msi_dcast$Bin_format, levels = c("Investigate Actionability","Cancer Gene Census","Cancer Hotspot","Cancer Pathway","Cosmic","Biologically Relevant"))
table(moalmanac_msi_dcast$Bin_format)
colnames(moalmanac_msi_dcast)[c(3,62)]<-c("Bin_Old","Bin")
# Upset plot Figure S8A
Figure_moalmanac_msi<-
  (
    upset(moalmanac_msi_dcast, callers, name=NULL, width_ratio=0.2, height_ratio = 0.8,  intersections = list(c("Union","MuTect2","Consensus2","VarScan2"),c("Union","Consensus2","Consensus3","MuTect2","MuSE","VarScan2","SomaticSniper"),c("Union","VarScan2"),c("Union","MuTect2"),c("Union","Consensus2","Consensus3","MuTect2","MuSE","VarScan2"),c("Union","Consensus2","Consensus3","MuSE","VarScan2","SomaticSniper"),c("Union","MuTect2","MuSE","Consensus2"),c("Union","MuSE"),c("Union","SomaticSniper")), keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("MC3"=(ggplot(mapping=aes(fill=MC3_Overlap))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "MC3 Overlap",labels=c("TRUE","FALSE",""),values=c('TRUE'="#2D708EFF", 'FALSE'="azure3"))+ylab("MC3 Overlap")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "VAF"=(ggplot(mapping=aes(y=tumor_f_mean))+geom_jitter(aes(color=total_coverage_mean_factor),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(1,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>130'="cadetblue1", '>75-130'="aquamarine1", '45-75'="antiquewhite3", '<45'="antiquewhite")))),
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),text=list(size=3.5),mapping=aes(fill=Bin)) + scale_fill_manual(values=c('Biologically Relevant'=cbPalette[1], 'Investigate Actionability'=cbPalette[2], 'Cancer Gene Census'=cbPalette[3],'Cancer Hotspot'=cbPalette[4],'Cancer Pathway'=cbPalette[5],'Cosmic'=cbPalette[6])) + theme(legend.key.size = unit(1,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Bin))) + geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=1500) +theme(legend.position = "none")+scale_fill_manual(values=c('Biologically Relevant'=cbPalette[1], 'Investigate Actionability'=cbPalette[2], 'Cancer Gene Census'=cbPalette[3],'Cancer Hotspot'=cbPalette[4],'Cancer Pathway'=cbPalette[5],'Cosmic'=cbPalette[6]))))
    + ggtitle("MSI Variants") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_moalmanac_msi, width = 9, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/FigureS8A.png",sep = "/"))
# Save results
saveRDS(Figure_moalmanac_msi,paste(dataDir,"Figures/FigureS8A.rds",sep = "/"))
# Write report
moalmanac_msi_dcast_report<-moalmanac_msi_dcast
moalmanac_msi_dcast_report[callers_moalmanac]<-ifelse(moalmanac_msi_dcast_report[callers_moalmanac]==TRUE,1,0)
writexl::write_xlsx(moalmanac_msi_dcast_report,paste(dataDir,"Supplementary_files/Supplementary_file_10_moalmanac_msi_dcast.xlsx",sep = "/"))

################################################################################
### Step 6: Plot Figure S8B
################################################################################

## Read supplementary data from papers: Bonneville et al. 2017 & Cortes-Ciriano et al. 2017
# Bonneville et al., 2017
Boneville<-as.data.frame(read_excel(paste(dataDir,"Resources/MOAlmanac_MSI/Bonneville_Suppl.xlsx",sep = "/")))
# Select TCGA cases
Boneville<-Boneville[grepl("TCGA",Boneville$`Case ID`),]
# Format data
Boneville$`Cancer Type`<-gsub("TCGA-","",Boneville$`Cancer Type`)
table(Boneville$`Cancer Type`)
# Define MSI-H status
Boneville$MSI_Status<-ifelse(Boneville$`MANTIS Score` > 0.4, "MSI-H",NA)
table(Boneville$`Cancer Type`,Boneville$MSI_Status)
# Cortes-Ciriano et al. 2017
Ciriano<-as.data.frame(read_excel(paste(dataDir,"Resources/MOAlmanac_MSI/Cortes_Ciriano_Suppl1.xlsx",sep = "/")))
# Format data
Ciriano$MSI_category_nb_from_TCGA_consortium<-toupper(Ciriano$MSI_category_nb_from_TCGA_consortium)
table(Ciriano$MSI_category_nb_from_TCGA_consortium)
head(Ciriano)
#Select UCEC STAD COAD READ cases 
MSI_cancers<-c("UCEC","STAD","COAD","READ")
Boneville<-Boneville[Boneville$`Cancer Type` %in% MSI_cancers,]
Ciriano<-Ciriano[Ciriano$Cancer_type %in% MSI_cancers,]
table(Boneville$`Cancer Type`)
table(Ciriano$Cancer_type)
# Select MSI-H cases
Boneville<-Boneville[Boneville$MSI_Status == "MSI-H",]
Ciriano<-Ciriano[Ciriano$MSI_category_nb_from_TCGA_consortium == "MSI-H",]

## Merge with MOAlmanac MSI results
# Read MOAlmanac results with purity/ploidy adjusted and aggregated by variant and sample and annotation data included
moalmanac_msi_final<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_results_adjusted_annotated.rds",sep = "/"))
head(moalmanac_msi_final)
# Remove patient id, tumor_f,total coverage,t_vaf_purity,t_vaf_purity_ploidy and t_vaf_purity_ploidy_final columns for Upset plot
colnames(moalmanac_msi_final)
moalmanac_msi_final<-moalmanac_msi_final[,-c(31,17:18,56:58)]
head(moalmanac_msi_final)
# Prepare data
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
Callers<-toupper(callers)
# Put caller column as last column
head(moalmanac_msi_final)
colnames(moalmanac_msi_final)
moalmanac_msi_final<-moalmanac_msi_final[,c(1:36,39:54,37:38)]
colnames(moalmanac_msi_final)[42:44]<-c("Cancer_DNA_fraction","Subclonal_genome_fraction","Genome_doublings")
# Transform blank spaces to NA!!!
moalmanac_msi_final[moalmanac_msi_final==""] <- NA
# Format patient id to merge
head(moalmanac_msi_final)
moalmanac_msi_final$tumor_sample_barcode_redux<-substr(moalmanac_msi_final$tumor_sample_barcode,1,12)
# Format column names to merge both files
head(moalmanac_msi_final)
colnames(Boneville)
colnames(Ciriano)
colnames(Boneville)[1]<-"tumor_sample_barcode_redux"
colnames(Ciriano)[2]<-"tumor_sample_barcode_redux"
# Select MSI-H samples
Boneville<-Boneville[!is.na(Boneville$tumor_sample_barcode_redux),]
Ciriano<-Ciriano[!is.na(Ciriano$tumor_sample_barcode_redux),]
sum(is.na(Ciriano$MSI_category_nb_from_TCGA_consortium))
length(unique(Boneville$tumor_sample_barcode_redux)) #348 MSI-H samples
length(unique(Ciriano$tumor_sample_barcode_redux)) #187 MSI-H samples
head(Boneville)
head(Ciriano)
MSI_H<-as.data.frame(unique(c(Boneville$tumor_sample_barcode_redux,Ciriano$tumor_sample_barcode_redux))) #355 unique MSI-H samples
MSI_H$MSI_Status<-"MSI-H"
colnames(MSI_H)[1]<-"tumor_sample_barcode_redux"
# Merge both files
moalmanac_msi_final_merge<-Reduce(function(x,y) merge(x,y,by="tumor_sample_barcode_redux",all=TRUE) ,list(moalmanac_msi_final, MSI_H))
dim(moalmanac_msi_final_merge)
head(moalmanac_msi_final_merge)
View(moalmanac_msi_final_merge)

## Plot figure S8B
# Prepare data for plot
msi_columns<-c("MSI_Status")
moalmanac_msi_final_merge[msi_columns][is.na(moalmanac_msi_final_merge[msi_columns])] <- 0
moalmanac_msi_final_merge[msi_columns][moalmanac_msi_final_merge[msi_columns] == "MSI-H"]<-1
View(moalmanac_msi_final_merge)
# Transform blank spaces to NA!!!
moalmanac_msi_final_merge[moalmanac_msi_final_merge==""] <- NA
# Put caller column as last column
head(moalmanac_msi_final_merge)
colnames(moalmanac_msi_final_merge)
moalmanac_msi_final_merge<-moalmanac_msi_final_merge[,c(1:53,56,54:55)]
# Save data
saveRDS(moalmanac_msi_final_merge,paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_final_merge.rds",sep = "/"))

## Plot Figure S8B
# Read data
moalmanac_msi_final_merge<-readRDS(paste(dataDir,"Resources/MOAlmanac_MSI/moalmanac_msi_final_merge.rds",sep = "/"))
#Select MSI-H tumors (UCEC,COAD,READ,STAD)
MSI_H_Cancers<-c("UCEC","COAD","READ","STAD")
moalmanac_msi_final_merge<-moalmanac_msi_final_merge[moalmanac_msi_final_merge$Cancer %in% MSI_H_Cancers,]
# Reshape data
moalmanac_msi_final_merge_dcast<-reshape2::dcast(moalmanac_msi_final_merge, tumor_sample_barcode_redux+Cancer+MSI_Status~Caller, fun.aggregate = length)
colnames(moalmanac_msi_final_merge_dcast)[c(3:10)]<-c("MSI-H","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers_msi<-c(callers,"MSI-H")
moalmanac_msi_final_merge_dcast[c(callers_msi)]<-moalmanac_msi_final_merge_dcast[c(callers_msi)] != 0
head(moalmanac_msi_final_merge_dcast)
# Upset plot Figure S8B
Figure_moalmanac_msi_benchmark<-
  (
    upset(moalmanac_msi_final_merge_dcast, callers_msi, name=NULL, width_ratio=0.25, height_ratio = 0.75,  intersections = list(c("Union","MuTect2","VarScan2","MSI-H","Consensus2"),c("Union","Consensus2","Consensus3","MuTect2","VarScan2","MuSE","SomaticSniper","MSI-H"),c("Union","Consensus2","Consensus3","MuTect2","VarScan2","MuSE","SomaticSniper"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuTect2"),c("Union","VarScan2","MSI-H"),c("Union","VarScan2"),c("Union","MuTect2","VarScan2","MSI-H"),c("Union","Consensus2","Consensus3","MuTect2","VarScan2","MuSE","MSI-H"),c("Union","MuTect2","MSI-H"),c("Union","MuSE"),c("Union","SomaticSniper")), keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),text=list(size=3),mapping=aes(fill=Cancer)) + scale_fill_manual(values=c('COAD'="darkgoldenrod1", 'READ'="darkgoldenrod4", 'STAD'="darkblue",'UCEC'="darkcyan")) + theme(legend.key.size = unit(1,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Cancer))) + geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.75)+ expand_limits(y=400) +theme(legend.position = "none")+scale_fill_manual(values=c('COAD'="darkgoldenrod1", 'READ'="darkgoldenrod4", 'STAD'="darkblue",'UCEC'="darkcyan"))))
    + ggtitle("Samples with MSI variants") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_moalmanac_msi_benchmark, width = 9, height = 6, dpi = 300, filename = paste(dataDir,"Figures/FigureS8B.png",sep = "/"))
# Save results
saveRDS(Figure_moalmanac_msi_benchmark,paste(dataDir,"Figures/FigureS8B.rds",sep = "/"))
# Write report
moalmanac_msi_final_merge_dcast_report<-moalmanac_msi_final_merge_dcast
moalmanac_msi_final_merge_dcast_report[callers_msi]<-ifelse(moalmanac_msi_final_merge_dcast_report[callers_msi]==TRUE,1,0)
writexl::write_xlsx(moalmanac_msi_final_merge_dcast_report,paste(dataDir,"Supplementary_files/Supplementary_file_X_moalmanac_msi_final_merge_dcast.xlsx",sep = "/"))







