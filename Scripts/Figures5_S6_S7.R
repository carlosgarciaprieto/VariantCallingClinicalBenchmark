#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figures 5, Figure S6 and Figure S7
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
### Step 1: Plot Figure S6
################################################################################

# Read MOAlmanac results with purity/ploidy adjusted and aggregated by variant and sample and annotation data included (see Prepare_data_Figures5_S6_S7.R)
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_results_adjusted_annotated.rds",sep = "/"))
# Remove patient id, tumor_f,total coverage,t_vaf_purity,t_vaf_purity_ploidy and t_vaf_purity_ploidy_final columns for Upset plot
colnames(moalmanac_final)
moalmanac_final<-moalmanac_final[,-c(12:13,41,61:63)]
head(moalmanac_final)
#Plot MOAlmanac results
callers<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
Callers<-toupper(callers)

## Prepare data
# Put caller column as last column
head(moalmanac_final)
colnames(moalmanac_final)
moalmanac_final<-moalmanac_final[,c(1:39,42:59,40:41)]
colnames(moalmanac_final)[47:49]<-c("Cancer_DNA_fraction","Subclonal_genome_fraction","Genome_doublings")
# Transform blank spaces to NA!!!
moalmanac_final[moalmanac_final==""] <- NA
# Reshape data
moalmanac_global_dcast<-reshape2::dcast(moalmanac_final, ...~Caller)
colnames(moalmanac_global_dcast)[c(4,59:65)]<-c("Bin","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
callers_moalmanac<-c("MC3_Overlap","GDC_Validation_Status",callers)
moalmanac_global_dcast$MC3_Overlap<-ifelse(moalmanac_global_dcast$MC3_Overlap == TRUE, 1, 0) #95.8% clinically actionable variants in MC3
moalmanac_global_dcast$GDC_Validation_Status<-ifelse(moalmanac_global_dcast$GDC_Validation_Status == "Valid", 1, 0) #8.5% clinically actionable variants have been validated
moalmanac_global_dcast[callers_moalmanac]<-moalmanac_global_dcast[callers_moalmanac] != 0
head(moalmanac_global_dcast)
moalmanac_global_dcast$MC3_Overlap<-factor(moalmanac_global_dcast$MC3_Overlap,levels = c("TRUE","FALSE"))
moalmanac_global_dcast$GDC_Validation_Status<-factor(moalmanac_global_dcast$GDC_Validation_Status,levels = c("TRUE","FALSE"))
#Save file
saveRDS(moalmanac_global_dcast,paste(dataDir,"/Resources/MOAlmanac/moalmanac_global_dcast.rds",sep = "/"))

## Plot all clinically actionable somatic variants (Figure S6)
# Set color palette
cbPalette <- c("#E69F00","#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Read data
moalmanac_global_dcast<-readRDS(paste(dataDir,"/Resources/MOAlmanac/moalmanac_global_dcast.rds",sep = "/"))
moalmanac_global_dcast_somatic<-moalmanac_global_dcast[moalmanac_global_dcast$feature_type=="Somatic Variant",]
head(moalmanac_global_dcast_somatic)
# Check coverage data
summary(moalmanac_global_dcast_somatic$total_coverage_mean) #Median Coverage 75 Mean Coverage 118
moalmanac_global_dcast_somatic$total_coverage_mean_factor<-ifelse(moalmanac_global_dcast_somatic$total_coverage_mean<45,"<45",
                                                                  ifelse(moalmanac_global_dcast_somatic$total_coverage_mean<=75,"45-75",
                                                                         ifelse(moalmanac_global_dcast_somatic$total_coverage_mean<=130,">75-130",">130")))
head(moalmanac_global_dcast_somatic)
moalmanac_global_dcast_somatic$total_coverage_mean_factor<-factor(moalmanac_global_dcast_somatic$total_coverage_mean_factor,levels = c(">130",">75-130","45-75","<45"))
table(moalmanac_global_dcast_somatic$total_coverage_mean_factor)
# Check stage information
moalmanac_global_dcast_somatic$ajcc_stage_parsed<-factor(moalmanac_global_dcast_somatic$ajcc_stage_parsed,levels = c("Stage I","Stage II","Stage III","Stage IV", "Not Available"))
table(moalmanac_global_dcast_somatic$ajcc_stage_parsed)
# Check categories
moalmanac_global_dcast_somatic$Bin<-factor(moalmanac_global_dcast_somatic$Bin, levels = c("Putatively Actionable","Investigate Actionability","Biologically Relevant"))
table(moalmanac_global_dcast_somatic$Bin)
# Upset plot
Figure_moalmanac_global<-
  (
    upset(moalmanac_global_dcast_somatic, callers, name=NULL, width_ratio=0.2, height_ratio = 0.8,  intersections = list(c("Union","MuTect2","VarScan2","MuSE","SomaticSniper","Consensus2","Consensus3"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuTect2"),c("Union","MuTect2","VarScan2","MuSE","Consensus2","Consensus3"),c("Union","VarScan2","MuSE","SomaticSniper","Consensus2","Consensus3"),c("Union","MuTect2","Consensus2","MuSE"),c("Union","VarScan2"),c("Union","MuTect2","VarScan2","SomaticSniper","Consensus2","Consensus3"),c("Union","MuSE"),c("Union","SomaticSniper")), keep_empty_groups = FALSE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("MC3"=(ggplot(mapping=aes(fill=MC3_Overlap))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "MC3 Overlap",labels=c("TRUE","FALSE",""),values=c('TRUE'="#2D708EFF", 'FALSE'="azure3"))+ylab("MC3 Overlap")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          "VAF"=(ggplot(mapping=aes(y=tumor_f_mean))+geom_jitter(aes(color=total_coverage_mean_factor),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(1,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>130'="cadetblue1", '>75-130'="aquamarine1", '45-75'="antiquewhite3", '<45'="antiquewhite")))),
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),text=list(size=3.5),mapping=aes(fill=Bin)) + scale_fill_manual(values=c('Biologically Relevant'=cbPalette[1], 'Investigate Actionability'=cbPalette[2], 'Putatively Actionable'=cbPalette[3])) + theme(legend.key.size = unit(1,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=Bin)))+geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=46000)+theme(legend.position = "none")+scale_fill_manual(values=c('Biologically Relevant'=cbPalette[1],'Investigate Actionability'=cbPalette[2],'Putatively Actionable'=cbPalette[3]))))
    + ggtitle("Clinically Actionable Somatic Variants") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_moalmanac_global, width = 9, height = 6.5, dpi = 300, filename = paste(dataDir,"Figures/FigureS6.png",sep = "/"))
# Save results
saveRDS(Figure_moalmanac_global,paste(dataDir,"Figures/Figure_moalmanac_global.rds",sep="/"))
# Write report
moalmanac_global_dcast_somatic_report<-moalmanac_global_dcast_somatic
moalmanac_global_dcast_somatic_report[callers_moalmanac]<-ifelse(moalmanac_global_dcast_somatic_report[callers_moalmanac]==TRUE,1,0)
writexl::write_xlsx(moalmanac_global_dcast_somatic_report,paste(dataDir,"Supplementary_files/Supplementary_file_6_moalmanac_global_dcast_somatic.xlsx",sep = "/"))

################################################################################
### Step 2: Plot Figure 5
################################################################################

# Read MOAlmanac results with purity/ploidy adjusted and aggregated by variant and sample and annotation data included (see Prepare_data_Figures5_S6_S7.R)
moalmanac_final<-readRDS(paste(dataDir,"Resources/MOAlmanac/moalmanac_results_adjusted_annotated.rds",sep = "/"))
# Remove patient id, tumor_f,total coverage,t_vaf_purity,t_vaf_purity_ploidy and t_vaf_purity_ploidy_final columns for Upset plot
colnames(moalmanac_final)
moalmanac_final<-moalmanac_final[,-c(12:13,41,61:63)]
head(moalmanac_final)
## Prepare data
# Put caller column as last column
head(moalmanac_final)
colnames(moalmanac_final)
moalmanac_final<-moalmanac_final[,c(1:39,42:59,40:41)]
colnames(moalmanac_final)[47:49]<-c("Cancer_DNA_fraction","Subclonal_genome_fraction","Genome_doublings")
# Transform blank spaces to NA!!!
moalmanac_final[moalmanac_final==""] <- NA
# Read reshape data
moalmanac_global_dcast<-readRDS(paste(dataDir,"/Resources/MOAlmanac/moalmanac_global_dcast.rds",sep = "/"))

## Prepare therapeutic sensitivity, resistance and prognostic data
moalmanac_final$favorable_prognosis<-ifelse(moalmanac_final$favorable_prognosis==1,"Favorable",ifelse(moalmanac_final$favorable_prognosis==0,"Unfavorable",NA))
moalmanac_global_dcast$favorable_prognosis<-ifelse(moalmanac_global_dcast$favorable_prognosis==1,"Favorable",ifelse(moalmanac_global_dcast$favorable_prognosis==0,"Unfavorable",NA))
colnames(moalmanac_global_dcast)
moalmanac_final$sensitive_predictive_implication<-factor(moalmanac_final$sensitive_predictive_implication,levels = c("FDA-Approved","Guideline","Clinical trial","Clinical evidence","Preclinical","Inferential"))
moalmanac_final$resistance_predictive_implication<-factor(moalmanac_final$resistance_predictive_implication,levels = c("FDA-Approved","Guideline","Clinical trial","Clinical evidence","Preclinical","Inferential"))
moalmanac_final$prognostic_predictive_implication<-factor(moalmanac_final$prognostic_predictive_implication,levels = c("FDA-Approved","Guideline","Clinical trial","Clinical evidence","Preclinical","Inferential"))

## Select those clinically actionable somatic variants in which the disease for which the association has been reported coincides with the cancer type of the tumor under analysis (including PANCANCER).
# Therapeutic sensitivity
moalmanac_sensitivity<-moalmanac_final[mapply(grepl, moalmanac_final$Cancer, moalmanac_final$sensitive_description_cancer) | mapply(grepl, "PANCANCER", moalmanac_final$sensitive_description_cancer),]
dim(moalmanac_sensitivity)
table(moalmanac_sensitivity$Cancer,moalmanac_sensitivity$sensitive_description_cancer)
# Therapeutic resistance
moalmanac_resistance<-moalmanac_final[mapply(grepl, moalmanac_final$Cancer, moalmanac_final$resistance_description_cancer) | mapply(grepl, "PANCANCER", moalmanac_final$resistance_description_cancer),]
dim(moalmanac_resistance)
table(moalmanac_resistance$Cancer,moalmanac_resistance$resistance_description_cancer)
# Disease prognosis
moalmanac_prognostic<-moalmanac_final[mapply(grepl, moalmanac_final$Cancer, moalmanac_final$prognostic_description_cancer) | mapply(grepl, "PANCANCER", moalmanac_final$prognostic_description_cancer),]
dim(moalmanac_prognostic)
table(moalmanac_prognostic$Cancer,moalmanac_prognostic$prognostic_description_cancer)

## Prepare data for plot
colnames(moalmanac_sensitivity)
moalmanac_sensitivity_dcast<-reshape2::dcast(moalmanac_sensitivity[moalmanac_sensitivity$feature_type=="Somatic Variant",], ... ~ Caller)
colnames(moalmanac_sensitivity_dcast)[c(5,59:65)]<-c("Evidence","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
moalmanac_resistance_dcast<-reshape2::dcast(moalmanac_resistance[moalmanac_resistance$feature_type=="Somatic Variant",], ... ~ Caller)
colnames(moalmanac_resistance_dcast)[c(6,59:65)]<-c("Evidence","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
moalmanac_prognostic_dcast<-reshape2::dcast(moalmanac_prognostic[moalmanac_prognostic$feature_type=="Somatic Variant",], ... ~ Caller)
colnames(moalmanac_prognostic_dcast)[c(7,59:65)]<-c("Evidence","Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
# Therapeutic sensitivity
callers_moalmanac<-c("MC3_Overlap","GDC_Validation_Status",callers)
moalmanac_sensitivity_dcast$MC3_Overlap<-ifelse(moalmanac_sensitivity_dcast$MC3_Overlap == TRUE, 1, 0) #92.5% clinically actionable sensitive variants in MC3
moalmanac_sensitivity_dcast$GDC_Validation_Status<-ifelse(moalmanac_sensitivity_dcast$GDC_Validation_Status == "Valid", 1, 0) #7.7% clinically actionable sensitive variants have been validated
moalmanac_sensitivity_dcast[callers_moalmanac]<-moalmanac_sensitivity_dcast[callers_moalmanac] != 0
head(moalmanac_sensitivity_dcast)
moalmanac_sensitivity_dcast$MC3_Overlap<-factor(moalmanac_sensitivity_dcast$MC3_Overlap,levels = c("TRUE","FALSE"))
moalmanac_sensitivity_dcast$GDC_Validation_Status<-factor(moalmanac_sensitivity_dcast$GDC_Validation_Status,levels = c("TRUE","FALSE"))
# Therapeutic resistance
moalmanac_resistance_dcast$MC3_Overlap<-ifelse(moalmanac_resistance_dcast$MC3_Overlap == TRUE, 1, 0) #88.5% clinically actionable resistance variants in MC3
moalmanac_resistance_dcast$GDC_Validation_Status<-ifelse(moalmanac_resistance_dcast$GDC_Validation_Status == "Valid", 1, 0) #7.4% clinically actionable resistance variants have been validated
moalmanac_resistance_dcast[callers_moalmanac]<-moalmanac_resistance_dcast[callers_moalmanac] != 0
head(moalmanac_resistance_dcast)
moalmanac_resistance_dcast$MC3_Overlap<-factor(moalmanac_resistance_dcast$MC3_Overlap,levels = c("TRUE","FALSE"))
moalmanac_resistance_dcast$GDC_Validation_Status<-factor(moalmanac_resistance_dcast$GDC_Validation_Status,levels = c("TRUE","FALSE"))
# Disease prognosis
moalmanac_prognostic_dcast$MC3_Overlap<-ifelse(moalmanac_prognostic_dcast$MC3_Overlap == TRUE, 1, 0) #91.8% clinically actionable prognostic variants in MC3
moalmanac_prognostic_dcast$GDC_Validation_Status<-ifelse(moalmanac_prognostic_dcast$GDC_Validation_Status == "Valid", 1, 0) #21.3% clinically actionable prognostic variants have been validated
moalmanac_prognostic_dcast[callers_moalmanac]<-moalmanac_prognostic_dcast[callers_moalmanac] != 0
head(moalmanac_prognostic_dcast)
moalmanac_prognostic_dcast$MC3_Overlap<-factor(moalmanac_prognostic_dcast$MC3_Overlap,levels = c("TRUE","FALSE"))
moalmanac_prognostic_dcast$GDC_Validation_Status<-factor(moalmanac_prognostic_dcast$GDC_Validation_Status,levels = c("TRUE","FALSE"))
# Save files
saveRDS(moalmanac_sensitivity_dcast,paste(dataDir,"/Resources/MOAlmanac/moalmanac_sensitivity_dcast.rds",sep = "/"))
saveRDS(moalmanac_resistance_dcast,paste(dataDir,"/Resources/MOAlmanac/moalmanac_resistance_dcast.rds",sep = "/"))
saveRDS(moalmanac_prognostic_dcast,paste(dataDir,"/Resources/MOAlmanac/moalmanac_prognostic_dcast.rds",sep = "/"))

## Prepare coverage data for Figure 5
# Therapeutic sensitivity
summary(moalmanac_sensitivity_dcast$total_coverage_mean) #Median Coverage 87 Mean Coverage 121
moalmanac_sensitivity_dcast$total_coverage_mean_factor<-ifelse(moalmanac_sensitivity_dcast$total_coverage_mean<50,"<50",
                                                               ifelse(moalmanac_sensitivity_dcast$total_coverage_mean<=90,"50-90",
                                                                      ifelse(moalmanac_sensitivity_dcast$total_coverage_mean<=150,">90-150",">150")))
head(moalmanac_sensitivity_dcast)
moalmanac_sensitivity_dcast$total_coverage_mean_factor<-factor(moalmanac_sensitivity_dcast$total_coverage_mean_factor,levels = c(">150",">90-150","50-90","<50"))
table(moalmanac_sensitivity_dcast$total_coverage_mean_factor)
# Therapeutic resistance
summary(moalmanac_resistance_dcast$total_coverage_mean) #Median Coverage 90 Mean Coverage 116
moalmanac_resistance_dcast$total_coverage_mean_factor<-ifelse(moalmanac_resistance_dcast$total_coverage_mean<50,"<50",
                                                              ifelse(moalmanac_resistance_dcast$total_coverage_mean<=90,"50-90",
                                                                     ifelse(moalmanac_resistance_dcast$total_coverage_mean<=150,">90-150",">150")))
head(moalmanac_resistance_dcast)
moalmanac_resistance_dcast$total_coverage_mean_factor<-factor(moalmanac_resistance_dcast$total_coverage_mean_factor,levels = c(">150",">90-150","50-90","<50"))
table(moalmanac_resistance_dcast$total_coverage_mean_factor)
# Disease prognosis
summary(moalmanac_prognostic_dcast$total_coverage_mean) #Median Coverage 96 Mean Coverage 133
moalmanac_prognostic_dcast$total_coverage_mean_factor<-ifelse(moalmanac_prognostic_dcast$total_coverage_mean<50,"<50",
                                                              ifelse(moalmanac_prognostic_dcast$total_coverage_mean<=90,"50-90",
                                                                     ifelse(moalmanac_prognostic_dcast$total_coverage_mean<=150,">90-150",">150")))
head(moalmanac_prognostic_dcast)
moalmanac_prognostic_dcast$total_coverage_mean_factor<-factor(moalmanac_prognostic_dcast$total_coverage_mean_factor,levels = c(">150",">90-150","50-90","<50"))
table(moalmanac_prognostic_dcast$total_coverage_mean_factor)

## Plot figure 5
# Plot first with n_intersections=10 and then plot specific intersections of interest i.e intersections=list()
Figure_moalmanac_clinical_sensitive_resistance<-
  (
    upset(moalmanac_sensitivity_dcast, callers, name=NULL, width_ratio=0.25, height_ratio = 0.75, intersections=list(c("Union","MuTect2","VarScan2","Consensus2","MuSE","Consensus3","SomaticSniper"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuTect2"),c("Union","MuTect2","VarScan2","Consensus2","MuSE","Consensus3"),c("Union","VarScan2","Consensus2","MuSE","Consensus3","SomaticSniper"),c("Union","VarScan2"),c("Union","Consensus2","MuTect2","MuSE"),c("Union","Consensus2","Consensus3","MuTect2","VarScan2","SomaticSniper"),c("Union","MuSE"),c("Union","SomaticSniper")),keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))), annotations=list("MC3"=(ggplot(mapping=aes(fill=MC3_Overlap))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "MC3 Overlap",labels=c("TRUE","FALSE",""),values=c('TRUE'="#2D708EFF", 'FALSE'="azure3"))+ylab("MC3 Overlap")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "VAF"=(ggplot(mapping=aes(y=tumor_f_mean))+geom_jitter(aes(color=total_coverage_mean_factor),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(1,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>150'="cadetblue1", '>90-150'="aquamarine1", '50-90'="antiquewhite3", '<50'="antiquewhite")))),
          
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=`Evidence`)) + scale_fill_manual(values=c('FDA-Approved'=cbPalette[1], 'Guideline'=cbPalette[2],'Clinical trial'=cbPalette[3], 'Clinical evidence'=cbPalette[4], 'Preclinical'=cbPalette[5], 'Inferential'=cbPalette[6])) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=`Evidence`))) + geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=7000) +theme(legend.position = "none")+scale_fill_manual(values=c(cbPalette[1],cbPalette[2],cbPalette[3],cbPalette[4],cbPalette[5],cbPalette[6]))))
    
    + ggtitle("Therapeutic sensitivity") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    +
      upset(moalmanac_resistance_dcast, callers, name=NULL, width_ratio=0.25, height_ratio = 0.75,  intersections=list(c("Union","MuTect2","VarScan2","Consensus2","MuSE","Consensus3","SomaticSniper"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuTect2"),c("Union","MuTect2","VarScan2","Consensus2","MuSE","Consensus3"),c("Union","VarScan2"),c("Union","VarScan2","SomaticSniper","Consensus2"),c("Union","Consensus2","MuTect2","MuSE"),c("Union","MuSE"),c("Union","Consensus2","Consensus3","MuTect2","VarScan2","SomaticSniper"),c("Union","SomaticSniper")), keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("MC3"=(ggplot(mapping=aes(fill=MC3_Overlap))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "MC3 Overlap",labels=c("TRUE","FALSE",""),values=c('TRUE'="#2D708EFF", 'FALSE'="azure3"))+ylab("MC3 Overlap")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "VAF"=(ggplot(mapping=aes(y=tumor_f_mean))+geom_jitter(aes(color=total_coverage_mean_factor),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(1,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>150'="cadetblue1", '>90-150'="aquamarine1", '50-90'="antiquewhite3", '<50'="antiquewhite")))), 
            
            base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=`Evidence`)) + scale_fill_manual(values=c('FDA-Approved'=cbPalette[1], 'Guideline'=cbPalette[2],'Clinical trial'=cbPalette[3], 'Clinical evidence'=cbPalette[4], 'Preclinical'=cbPalette[5], 'Inferential'=cbPalette[6])) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=`Evidence`)))+ geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=625)+theme(legend.position = "none")+scale_fill_manual(values=c(cbPalette[2],cbPalette[4],cbPalette[5]))))
    
    + ggtitle("Therapeutic resistance") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
ggsave(plot = Figure_moalmanac_clinical_sensitive_resistance, width = 11.5, height = 7, dpi = 300, filename = paste(dataDir,"Figures/Figure5.png",sep = "/"))
# Save results
saveRDS(Figure_moalmanac_clinical_sensitive_resistance,paste(dataDir,"Figures/Figure5.rds",sep = "/"))
# Write report
moalmanac_sensitivity_dcast_report<-moalmanac_sensitivity_dcast
moalmanac_sensitivity_dcast_report[callers_moalmanac]<-ifelse(moalmanac_sensitivity_dcast_report[callers_moalmanac]==TRUE,1,0)
writexl::write_xlsx(moalmanac_sensitivity_dcast_report,paste(dataDir,"Supplementary_files/Supplementary_file_7_moalmanac_sensitivity_dcast.xlsx",sep = "/"))
moalmanac_resistance_dcast_report<-moalmanac_resistance_dcast
moalmanac_resistance_dcast_report[callers_moalmanac]<-ifelse(moalmanac_resistance_dcast_report[callers_moalmanac]==TRUE,1,0)
writexl::write_xlsx(moalmanac_resistance_dcast_report,paste(dataDir,"Supplementary_files/Supplementary_file_8_moalmanac_resistance_dcast.xlsx",sep = "/"))

################################################################################
### Step 3: Plot Figure S7
################################################################################

## Plot figure S7
Figure_moalmanac_clinical_prognostic<-
  (
    upset(moalmanac_prognostic_dcast, callers, name=NULL, width_ratio=0.25, height_ratio = 0.75,  intersections = list(c("Union","Consensus3","Consensus2","MuTect2","VarScan2","MuSE","SomaticSniper"),c("Union","MuTect2","VarScan2","Consensus2"),c("Union","MuTect2"),c("Union","Consensus3","Consensus2","MuTect2","VarScan2","MuSE"),c("Union","Consensus2","MuTect2","MuSE"),c("Union","Consensus3","Consensus2","MuTect2","VarScan2","SomaticSniper"),c("Union","Consensus3","Consensus2","VarScan2","MuSE","SomaticSniper"),c("Union","VarScan2"),c("Union","MuSE"),c("Union","SomaticSniper")), keep_empty_groups = TRUE, wrap = TRUE, themes = upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=14,face = "bold")),'overall_sizes'=theme(axis.text.x=element_text(angle=90)))),annotations=list("MC3"=(ggplot(mapping=aes(fill=MC3_Overlap))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "MC3 Overlap",labels=c("TRUE","FALSE",""),values=c('TRUE'="#2D708EFF", 'FALSE'="azure3"))+ylab("MC3 Overlap")),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "VAF"=(ggplot(mapping=aes(y=tumor_f_mean))+geom_jitter(aes(color=total_coverage_mean_factor),size=0.25,na.rm=TRUE)+geom_violin(alpha=0.5,na.rm=TRUE)+theme(legend.key.size = unit(1,"line"))+guides(color = guide_legend(override.aes = list(size = 3)))+scale_color_manual(name="Coverage",values=c('>150'="cadetblue1", '>90-150'="aquamarine1", '50-90'="antiquewhite3", '<50'="antiquewhite"))),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       "Prognosis"=(ggplot(mapping=aes(fill=favorable_prognosis))+geom_bar(stat='count', position='fill',na.rm = TRUE)+theme(legend.key.size = unit(1,"line"))+scale_y_continuous(labels=scales::percent_format())+scale_fill_manual(name = "Prognosis",labels=c("Favorable","Unfavorable"),values=c('Favorable'="khaki1", 'Unfavorable'="khaki4"))+ylab("Prognosis"))), 
          
          base_annotations=list('Intersection size'=intersection_size(counts=TRUE,text_colors=c(on_bar="black",on_background="black"),mapping=aes(fill=`Evidence`)) + scale_fill_manual(values=c('FDA-Approved'=cbPalette[1], 'Guideline'=cbPalette[2],'Clinical trial'=cbPalette[3], 'Clinical evidence'=cbPalette[4], 'Preclinical'=cbPalette[5], 'Inferential'=cbPalette[6])) + theme(legend.key.size = unit(.85,"line"))), guides = "over", set_sizes = (upset_set_size(geom=geom_bar(aes(fill=`Evidence`)))+ geom_text(aes(label=..count..), hjust=1.1, stat='count', size = 2.5)+ expand_limits(y=1000) +theme(legend.position = "none")+scale_fill_manual(values=c(cbPalette[2],cbPalette[3],cbPalette[4],cbPalette[6]))))
    + ggtitle("Disease prognosis") + theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(plot = Figure_moalmanac_clinical_prognostic, width = 8, height = 8, dpi = 300, filename = paste(dataDir,"Figures/FigureS7.png",sep = "/"))
# Save results
saveRDS(Figure_moalmanac_clinical_prognostic,paste(dataDir,"Figures/FigureS7.rds",sep = "/"))
# Write report
moalmanac_prognostic_dcast_report<-moalmanac_prognostic_dcast
moalmanac_prognostic_dcast_report[callers_moalmanac]<-ifelse(moalmanac_prognostic_dcast_report[callers_moalmanac]==TRUE,1,0)
writexl::write_xlsx(moalmanac_prognostic_dcast_report,paste(dataDir,"Supplementary_files/Supplementary_file_9_moalmanac_prognostic_dcast.xlsx",sep = "/"))
#===============================================================================