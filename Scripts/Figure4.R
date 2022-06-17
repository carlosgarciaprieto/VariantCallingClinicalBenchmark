#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Figure 4
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
### Step 1: Run deconstructSigs
################################################################################

## Set the list of 5 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
cancer<-c("ACC","BLCA","BRCA","PRAD","UCEC")
caller<-c("CONSENSUS2","CONSENSUS3","MUSE","MUTECT2","SOMATICSNIPER","UNION","VARSCAN2")
Caller<-c("Consensus2","Consensus3","MuSE","MuTect2","SomaticSniper","Union","VarScan2")
  
# Load COSMIC reference matrix
load(paste(dataDir,"MutationalSignatures/signatures.genome.cosmic.v3.may2019.rda",sep = "/"))

# Select cancer specific mutational signatures expected to be found according to literature (Alexandrov et al.)
ACC_signatures<-c("SBS1","SBS2","SBS5","SBS13","SBS29","SBS39","SBS40")
BLCA_signatures<-c("SBS1","SBS2","SBS5","SBS8","SBS13","SBS29","SBS40")
BRCA_signatures<-c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS9","SBS13","SBS17a","SBS17b","SBS18","SBS37","SBS40","SBS41")
PRAD_signatures<-c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS13","SBS18","SBS33","SBS37","SBS40","SBS41")
UCEC_signatures<-c("SBS1","SBS2","SBS3","SBS5","SBS6","SBS9","SBS10a","SBS10b","SBS13","SBS14","SBS15","SBS18","SBS26","SBS28","SBS40","SBS41")
signatures_list<-list(ACC_signatures,BLCA_signatures,BRCA_signatures,PRAD_signatures,UCEC_signatures)

# For each of the 5 TCGA cancer types, we downloaded the four different Somatic aggregated MAF files with all the somatic mutations for each variant caller (MuSE, MuTect2, SomaticSniper and VarScan2), created the Consensus of 2, Consensus of 3 and Union files, and saved them under Maf_files folder
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  for(k in 1:length(caller)){
    iteration<-NULL
    maf<-NULL
    maf_df<-NULL
    sigs.input<-NULL
    results<-NULL
    expo<-NULL
    Signature.unknown<-NULL
    o<-NULL
    iteration<-paste(cancer[i],caller[k],sep = "_")
    print(iteration)
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    maf_df<-as.data.frame(maf@data)
    sigs.input <- mut.to.sigs.input(mut.ref = maf_df, 
                                    sample.id = "Tumor_Sample_Barcode", 
                                    chr = "Chromosome", 
                                    pos = "Start_Position", 
                                    ref = "Reference_Allele", 
                                    alt = "Tumor_Seq_Allele2",
                                    bsg = BSgenome.Hsapiens.UCSC.hg38)
    # filter samples with few mutations
    sigs.input <- sigs.input[rowSums(sigs.input)>50,]
    # initialize a list of the length of samples 
    results <- vector("list", nrow(sigs.input))
    names(results) <- row.names(sigs.input)
    # run the estimation of exposures for each sample and save the results in the list
    for( sID in row.names(sigs.input) ){
      results[[sID]] <- whichSignatures(sigs.input, # the matrix generated with mut.to.sigs.input 
                                        sample.id=sID, # the current sample ID
                                        signatures.ref=signatures.genome.cosmic.v3.may2019, # the data.frame with the signatures that comes with deconstructSigs
                                        associated=signatures_list[[i]] , #signatures expected to be found
                                        tri.counts.method="exome2genome", # which normalization method to use
                                        contexts.needed=TRUE) # set to TRUE if your input matrix contains counts instead of frequencies
    }
    # convert the exposures for each sample into a sample x signatures matrix
    expo <- do.call("rbind", sapply(results, "[", 1))
    # add the unknown value to the matrix such that the contributions add up to 1 per sample
    Signature.unknown <- unlist(sapply(results, "[", 5))
    expo <- cbind(expo, Signature.unknown)
    # reorder samples by similarity in their signature profiles
    o <- row.names(expo)[hclust(dist(expo))$order]
    expo$Caller<-Caller[k]
    expo$Sample<-rownames(expo)
    expo$Sample<-gsub("weights","",expo$Sample)
    expo$Sample<-gsub("\\.","",expo$Sample)
    rownames(expo)<-NULL
    expo$Cancer<-cancer[i]
    saveRDS(expo, file = paste0(dataDir,"/MutationalSignatures/",Caller[k],"_",cancer[i],"_","signatures.rds"))
  }
}

################################################################################
### Step 2: Aggregate results by cancer type
################################################################################

for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  file_cancer<-NULL
  signatures<-NULL
  for(k in 1:length(caller)){
    file<-NULL
    file<-readRDS(paste0(dataDir,"/MutationalSignatures/",Caller[k],"_",cancer[i],"_","signatures.rds"))
    file_cancer<-rbind(file_cancer,file)
  }
  signatures<-reshape2::melt(file_cancer, id.vars = c("Caller","Cancer","Sample"))
  colnames(signatures)[4:5]<-c("Signature","Weight")
  saveRDS(signatures,file = paste0(dataDir,"/MutationalSignatures/",cancer[i],"_signatures.rds"))
}

################################################################################
### Step 3: Compute proportion of samples with signature
################################################################################

for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  file_signatures<-NULL
  file_signatures_proportion<-NULL
  file_signatures_final<-NULL
  file_signatures<-readRDS(paste0(dataDir,"/MutationalSignatures/",cancer[i],"_","signatures.rds"))
  #Compute proportions of tumor with the signature
  file_signatures_proportion<-aggregate(file_signatures$Weight, list(file_signatures$Signature,file_signatures$Caller), FUN=function(x) length(which((x)>0.06))/length(x)) #Set 0.06 as minimun threshold (deconstructSigs default)
  colnames(file_signatures_proportion)<-c("Signature","Caller","Proportion")
  file_signatures_proportion$Signature_Caller<-paste(file_signatures_proportion$Signature,file_signatures_proportion$Caller,sep = "_")
  #Merge results
  file_signatures$Signature_Caller<-paste(file_signatures$Signature,file_signatures$Caller,sep = "_")
  file_signatures_final<-merge(file_signatures,file_signatures_proportion[,c(3:4)],by="Signature_Caller")
  saveRDS(file_signatures_final,file = paste0(dataDir,"/MutationalSignatures/",cancer[i],"_signatures_proportion_final.rds"))
}

################################################################################
### Step 4: Plot Figure 4
################################################################################

## Aggregate results
aggregate_signatures<-NULL
for(i in 1:length(cancer)){
  print(i)
  print(cancer[i])
  file_signatures_final<-NULL
  file_signatures_final<-readRDS(file = paste0(dataDir,"/MutationalSignatures/",cancer[i],"_signatures_proportion_final.rds"))
  aggregate_signatures<-rbind(aggregate_signatures,file_signatures_final)
}  
# Save results
saveRDS(aggregate_signatures,file = paste0(dataDir,"/MutationalSignatures/","aggregate_signatures.rds"))

## Compute median weight per cancer and variant calling strategy
aggregate_signatures$Signature_Caller_Cancer<-paste(aggregate_signatures$Signature_Caller,aggregate_signatures$Cancer,sep = "_")
head(aggregate_signatures)
Median_weight<-aggregate(aggregate_signatures$Weight, list(aggregate_signatures$Signature,aggregate_signatures$Caller,aggregate_signatures$Cancer), FUN=median) 
colnames(Median_weight)<-c("Signature","Caller","Cancer","Median_Weight")
head(Median_weight)
Median_weight$Signature_Caller_Cancer<-paste(Median_weight$Signature,Median_weight$Caller,Median_weight$Cancer,sep = "_")
# Merge information
aggregate_signatures_plot<-merge(aggregate_signatures,Median_weight[,c(4:5)],by="Signature_Caller_Cancer")
dim(aggregate_signatures_plot)
dim(aggregate_signatures_plot)
head(aggregate_signatures_plot)
# Save data
saveRDS(aggregate_signatures_plot,paste(dataDir,"MutationalSignatures/aggregate_signatures_plot.rds",sep = "/"))
# Write report
aggregate_signatures_plot<-aggregate_signatures_plot[order(aggregate_signatures_plot$Cancer,aggregate_signatures_plot$Caller,aggregate_signatures_plot$Signature),]
writexl::write_xlsx(aggregate_signatures_plot,paste(dataDir,"Supplementary_files/Supplementary_file_5_Mutational_signatures_results_median_weight.xlsx",sep = "/"))

## Plot Figure 4
# Prepare data for plot (remove unknown signature contribution)
aggregate_signatures<-aggregate_signatures[aggregate_signatures$Signature!="Signature.unknown",]
aggregate_signatures$Signature<-factor(aggregate_signatures$Signature)
# Upset plot
mutsign_bubble<-ggplot(aggregate_signatures, aes(x = Caller, y = Signature)) +theme_minimal(base_size = 12)+ geom_tile(colour= "gray20", fill= "white", size = 0.1,stat = "identity" ) +
  geom_point(aes(size = Proportion, fill = Median_Weight), shape = 21)+ facet_grid (~as.factor (aggregate_signatures$Cancer), scales = "free_x") + 
  ylab (element_blank()) + xlab (element_blank()) + scale_size_continuous(range = c(-0.4,3.1),breaks = c(0,0.25,0.5,0.75,1)) +
  theme(axis.text.x = element_text (angle = 90, hjust = 1, vjust = 0.5),legend.key.size = unit(1,"line"),legend.text = element_text(size = 8, colour = "black", face = "bold")) +
  theme(strip.background = element_rect(colour="black", fill="antiquewhite",size=0.5, linetype="solid"),strip.text.x = element_text(size=10, color="black",face="bold")) + theme(axis.text.y = element_text(color="black", size=7, face="bold"),axis.text.x = element_text(color="black", size=7.5, face="bold")) +
  scale_y_discrete(limits = rev(levels(aggregate_signatures$Signature)))+ scale_fill_viridis(na.value = "white",option="A", begin = 0, end = 1, direction = -1) 

ggsave(plot = mutsign_bubble, width = 7, height = 7, dpi = 300,filename = paste(dataDir,"Figures/Figure4.png",sep = "/"))
# Save plot data
saveRDS(mutsign_bubble,paste(dataDir,"Figures/Figure4.rds",sep = "/"))
#===============================================================================