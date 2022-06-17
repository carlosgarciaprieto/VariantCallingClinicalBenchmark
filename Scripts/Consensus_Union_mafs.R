#!/usr/bin/env Rscript
################################################################################
# Detection of oncogenic and clinically actionable mutations in cancer genomes critically depends on variant calling tools
# Author: Carlos A. Garcia-Prieto
# Description: script used to generate Consensus2, Consensus3 and Union MAF files
# Usage: provide cancer type to produce derived MAF files as first argument and data directory as second argument (TCGA MAF files should be located under Maf_files folder in data directory)
# Usage: from the command line run Rscript --vanilla Consensus_Union_mafs.R Cancer_type Directory
################################################################################
### Define arguments 
args = commandArgs(trailingOnly=TRUE)
# Test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("One TCGA cancer type (first argument) and data directory (second argument) containing Maf_files folder with TCGA downloaded MAF files must be supplied", call.=FALSE)
} 

### Libraries
print("Loading libraries")
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
## Set cancer type and working directory
cancer<-args[1]
dataDir<-args[2]
setwd(dataDir)
#===============================================================================
### MAIN.
################################################################################
### Step 1: Read TCGA MAF files downloaded from GDC Data Portal
################################################################################
print(paste("¡¡¡Starting analysis of", cancer, "!!!", sep = " "))
## Set the list of 33 TCGA projects (cancer types) and 7 variant calling strategies to be analysed
#cancer<-c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
          #"ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML",
          #"LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
          #"PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT",
          #"THCA","THYM","UCEC","UCS","UVM")

caller<-c("MUSE","MUTECT2","SOMATICSNIPER","VARSCAN2")

## Read MAF files and write file with variants detected by each variant caller
print("Read MAF files and write file with variants detected by each variant caller")
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    maf_df<-NULL
    maf_df_final<-NULL
    filename<-NULL
    #Set filename
    filename<-paste(cancer[i],caller[k],sep=".")
    #Read maf
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    #Create unique variant ID
    maf@data$name<-paste(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode, sep = ":")
    #Format data
    maf@data$caller<-1
    colnames(maf@data)[122]<-paste0(caller[k])
    #Create data.frame
    maf_df<-data.frame(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode)
    maf_df$name<- paste(maf_df$maf.data.Chromosome, maf_df$maf.data.Start_Position, maf_df$maf.data.End_Position, maf_df$maf.data.Tumor_Seq_Allele1, maf_df$maf.data.Tumor_Seq_Allele2, maf_df$maf.data.Tumor_Sample_Barcode, sep = ":")
    maf_df$Caller<-1
    colnames(maf_df)[8]<-paste0(caller[k])
    maf_df_final<-maf_df[,c("name",paste0(caller[k]))]
    write.table(maf_df_final, paste0(dataDir,"/Maf_files/",filename,".txt"), row.names = F, col.names = T, quote = F,sep="\t")
  }
}

################################################################################
### Step 2: Compare variants detected by each variant caller per cancer type
################################################################################
print("Compare variants detected by each variant caller per cancer type")
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  MUSE<-NULL
  MUTECT2<-NULL
  SOMATICSNIPER<-NULL
  VARSCAN2<-NULL
  comparison<-NULL
  comparison_morethan1callers<-NULL
  comparison_morethan2callers<-NULL
  #Read variants detected by each caller
  MUSE<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUSE.txt"))
  MUTECT2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUTECT2.txt"))
  SOMATICSNIPER<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".SOMATICSNIPER.txt"))
  VARSCAN2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".VARSCAN2.txt"))
  #Compare variants detected by each variant caller and select those variants detected with more than one and two callers
  comparison<-Reduce(function(x,y) merge(x,y,by="name",all=TRUE) ,list(MUSE,MUTECT2,SOMATICSNIPER,VARSCAN2))  
  comparison[is.na(comparison)] <- 0
  comparison$sum<-apply(comparison[,2:5],1,function(x) sum(x))
  comparison_morethan1callers<-comparison[comparison$sum > 1,]
  comparison_morethan2callers<-comparison[comparison$sum > 2,]
  #Write comparison results
  write.table(comparison, paste0(dataDir,"/Maf_files/",cancer[i],"_comparison.txt"), row.names = F, col.names = T, quote = F,sep="\t")
  write.table(comparison_morethan1callers, paste0(dataDir,"/Maf_files/",cancer[i],"_comparison_morethan1callers.txt"), row.names = F, col.names = T, quote = F,sep="\t")
  write.table(comparison_morethan2callers, paste0(dataDir,"/Maf_files/",cancer[i],"_comparison_morethan2callers.txt"), row.names = F, col.names = T, quote = F,sep="\t")
}

################################################################################
### Step 3: Create Consensus of at least three variant callers files (CONSENSUS3)
################################################################################
print("Create Consensus of at least three variant callers files (CONSENSUS3)")
## Read variants detected by at least three variant callers 
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    maf_morethan2<-NULL
    filename<-NULL
    morethan2<-NULL
    #Set filename
    filename<-paste(cancer[i],caller[k],sep=".")
    #Read maf
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    #Create unique variant ID
    maf@data$name<-paste(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode, sep = ":")
    #Format data
    maf@data$caller<-1
    colnames(maf@data)[122]<-paste0(caller[k])
    #Select Consensus of at least three variant callers
    morethan2<-read.table(paste0(dataDir,"/Maf_files/",cancer[i],"_comparison_morethan2callers.txt"), header = T, sep = "\t")
    maf_morethan2<-merge(maf@data, comparison_morethan2callers[,c(1,6)], by = "name")
    #Format data
    maf_morethan2$src_vcf_id<-"NA"
    maf_morethan2$Caller<-paste0(caller[k])
    maf_morethan2<-maf_morethan2[,-122]
    #Write results
    write.table(maf_morethan2, paste0(dataDir,"/Maf_files/",filename,".morethan2.txt"), row.names = F, col.names = T, quote = F,sep="\t")
  }
}

## Prepare CONSENSUS3 MAF file
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  MUSE_morethan2<-NULL
  MUTECT2_morethan2<-NULL
  SOMATICSNIPER_morethan2<-NULL
  VARSCAN2_morethan2<-NULL
  Total_morethan2<-NULL
  Total.dt<-NULL
  CONSENSUS3<-NULL
  #Read variants detected by each caller that are also detected by at least two other variant callers
  MUSE_morethan2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUSE.morethan2.txt"))
  MUTECT2_morethan2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUTECT2.morethan2.txt"))
  SOMATICSNIPER_morethan2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".SOMATICSNIPER.morethan2.txt"))
  VARSCAN2_morethan2<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".VARSCAN2.morethan2.txt"))
  #Combine all variants detected by at least three variant callers
  Total_morethan2<-rbind(MUSE_morethan2,MUTECT2_morethan2,SOMATICSNIPER_morethan2,VARSCAN2_morethan2)
  #Prepare CONSENSUS3 MAF file
  Total.dt <- data.table(Total_morethan2)
  setkey(Total.dt, name)
  Total.dt$t_depth<-as.numeric(as.character(Total.dt$t_depth))
  Total.dt$t_ref_count<-as.numeric(as.character(Total.dt$t_ref_count))
  Total.dt$t_alt_count<-as.numeric(as.character(Total.dt$t_alt_count))
  Total.dt$n_depth<-as.numeric(as.character(Total.dt$n_depth))
  #Define mean reference and variant alleles read depths
  CONSENSUS3<-data.frame(Total.dt[, c(list(Hugo_Symbol=Hugo_Symbol[1],Entrez_Gene_Id=Entrez_Gene_Id[1],Center=Center[1],NCBI_Build=NCBI_Build[1],Chromosome=Chromosome[1],Start_Position=Start_Position[1],End_Position=End_Position[1],Strand=Strand[1],Variant_Classification=Variant_Classification[1],Variant_Type=Variant_Type[1],Reference_Allele=Reference_Allele[1],Tumor_Seq_Allele1=Tumor_Seq_Allele1[1],Tumor_Seq_Allele2=Tumor_Seq_Allele2[1],dbSNP_RS=dbSNP_RS[1],dbSNP_Val_Status=dbSNP_Val_Status[1],Tumor_Sample_Barcode=Tumor_Sample_Barcode[1],Matched_Norm_Sample_Barcode=Matched_Norm_Sample_Barcode[1],Match_Norm_Seq_Allele1=Match_Norm_Seq_Allele1[1],Match_Norm_Seq_Allele2=Match_Norm_Seq_Allele2[1],Tumor_Validation_Allele1=Tumor_Validation_Allele1[1],Tumor_Validation_Allele2=Tumor_Validation_Allele2[1],Match_Norm_Validation_Allele1=Match_Norm_Validation_Allele1[1],Match_Norm_Validation_Allele2=Match_Norm_Validation_Allele2[1],Verification_Status=Verification_Status[1],Validation_Status=Validation_Status[1],Mutation_Status=Mutation_Status[1],Sequencing_Phase=Sequencing_Phase[1],Sequence_Source=Sequence_Source[1],Validation_Method=Validation_Method[1],Score=Score[1],BAM_File=BAM_File[1],Sequencer=Sequencer[1],Tumor_Sample_UUID=Tumor_Sample_UUID[1],Matched_Norm_Sample_UUID=Matched_Norm_Sample_UUID[1],HGVSc=HGVSc[1],HGVSp=HGVSp[1],HGVSp_Short=HGVSp_Short[1],Transcript_ID=Transcript_ID[1],Exon_Number=Exon_Number[1],n_ref_count=n_ref_count[1],n_alt_count=n_alt_count[1],all_effects=all_effects[1],Allele=Allele[1],Gene=Gene[1],Feature=Feature[1],Feature_type=Feature_type[1],One_Consequence=One_Consequence[1],Consequence=Consequence[1],cDNA_position=cDNA_position[1],CDS_position=CDS_position[1],Protein_position=Protein_position[1],Amino_acids=Amino_acids[1],Codons=Codons[1],Existing_variation=Existing_variation[1],ALLELE_NUM=ALLELE_NUM[1],DISTANCE=DISTANCE[1],TRANSCRIPT_STRAND=TRANSCRIPT_STRAND[1],SYMBOL=SYMBOL[1],SYMBOL_SOURCE=SYMBOL_SOURCE[1],HGNC_ID=HGNC_ID[1],BIOTYPE=BIOTYPE[1],CANONICAL=CANONICAL[1],CCDS=CCDS[1],ENSP=ENSP[1],SWISSPROT=SWISSPROT[1],TREMBL=TREMBL[1],UNIPARC=UNIPARC[1],RefSeq=RefSeq[1],SIFT=SIFT[1],PolyPhen=PolyPhen[1],EXON=EXON[1],INTRON=INTRON[1],DOMAINS=DOMAINS[1],GMAF=GMAF[1],AFR_MAF=AFR_MAF[1],AMR_MAF=AMR_MAF[1],ASN_MAF=ASN_MAF[1],EAS_MAF=EAS_MAF[1],EUR_MAF=EUR_MAF[1],SAS_MAF=SAS_MAF[1],AA_MAF=AA_MAF[1],EA_MAF=EA_MAF[1],CLIN_SIG=CLIN_SIG[1],SOMATIC=SOMATIC[1],PUBMED=PUBMED[1],MOTIF_NAME=MOTIF_NAME[1],MOTIF_POS=MOTIF_POS[1],HIGH_INF_POS=HIGH_INF_POS[1],MOTIF_SCORE_CHANGE=MOTIF_SCORE_CHANGE[1],IMPACT=IMPACT[1],PICK=PICK[1],VARIANT_CLASS=VARIANT_CLASS[1],TSL=TSL[1],HGVS_OFFSET=HGVS_OFFSET[1],PHENO=PHENO[1],MINIMISED=MINIMISED[1],ExAC_AF=ExAC_AF[1],ExAC_AF_Adj=ExAC_AF_Adj[1],ExAC_AF_AFR=ExAC_AF_AFR[1],ExAC_AF_AMR=ExAC_AF_AMR[1],ExAC_AF_EAS=ExAC_AF_EAS[1],ExAC_AF_FIN=ExAC_AF_FIN[1],ExAC_AF_NFE=ExAC_AF_NFE[1],ExAC_AF_OTH=ExAC_AF_OTH[1],ExAC_AF_SAS=ExAC_AF_SAS[1],GENE_PHENO=GENE_PHENO[1],FILTER=FILTER[1],CONTEXT=CONTEXT[1],src_vcf_id=src_vcf_id[1],tumor_bam_uuid=tumor_bam_uuid[1],normal_bam_uuid=normal_bam_uuid[1],case_id=case_id[1],GDC_FILTER=GDC_FILTER[1],COSMIC=COSMIC[1],MC3_Overlap=MC3_Overlap[1],GDC_Validation_Status=GDC_Validation_Status[1],sum=sum[1]), lapply(.SD,mean)), by ="name", .SDcols=c("t_depth","t_ref_count","t_alt_count","n_depth")])
  CONSENSUS3$t_depth<-round(CONSENSUS3$t_depth)
  CONSENSUS3$t_ref_count<-round(CONSENSUS3$t_ref_count)
  CONSENSUS3$t_alt_count<-round(CONSENSUS3$t_alt_count)
  CONSENSUS3$n_depth<-round(CONSENSUS3$n_depth)
  #Final CONSENSUS3 MAF file
  CONSENSUS3<-CONSENSUS3[,c(1:40,119:122,41:118)]
  CONSENSUS3<-CONSENSUS3[,-c(1,122)]
  CONSENSUS3<-CONSENSUS3[order(CONSENSUS3$Chromosome, CONSENSUS3$Start_Position), ]
  #Write CONSENSUS3 MAF file
  write.table(CONSENSUS3, paste0(dataDir,"/Maf_files/TCGA_",cancer[i],"_CONSENSUS3.maf"), row.names = F, col.names = T, quote = F,sep="\t")
}

################################################################################
### Step 4: Create Consensus of at least two variant callers files (CONSENSUS2)
################################################################################
print("Create Consensus of at least two variant callers files (CONSENSUS2)")
## Read variants detected by at least two variant callers 
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    maf_morethan1<-NULL
    filename<-NULL
    morethan1<-NULL
    #Set filename
    filename<-paste(cancer[i],caller[k],sep=".")
    #Read maf
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    #Create unique variant ID
    maf@data$name<-paste(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode, sep = ":")
    #Format data
    maf@data$caller<-1
    colnames(maf@data)[122]<-paste0(caller[k])
    #Select Consensus of at least two variant callers
    morethan1<-read.table(paste0(dataDir,"/Maf_files/",cancer[i],"_comparison_morethan1callers.txt"), header = T, sep = "\t")
    maf_morethan1<-merge(maf@data, comparison_morethan1callers[,c(1,6)], by = "name")
    #Format data
    maf_morethan1$src_vcf_id<-"NA"
    maf_morethan1$Caller<-paste0(caller[k])
    maf_morethan1<-maf_morethan1[,-122]
    #Write results
    write.table(maf_morethan1, paste0(dataDir,"/Maf_files/",filename,".morethan1.txt"), row.names = F, col.names = T, quote = F,sep="\t")
  }
}

## Prepare CONSENSUS2 MAF file
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  MUSE_morethan1<-NULL
  MUTECT2_morethan1<-NULL
  SOMATICSNIPER_morethan1<-NULL
  VARSCAN2_morethan1<-NULL
  Total_morethan1<-NULL
  Total.dt<-NULL
  CONSENSUS2<-NULL
  #Read variants detected by each caller that are also detected by at least another variant caller
  MUSE_morethan1<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUSE.morethan1.txt"))
  MUTECT2_morethan1<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUTECT2.morethan1.txt"))
  SOMATICSNIPER_morethan1<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".SOMATICSNIPER.morethan1.txt"))
  VARSCAN2_morethan1<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".VARSCAN2.morethan1.txt"))
  #Combine all variants detected by at least three variant callers
  Total_morethan1<-rbind(MUSE_morethan1,MUTECT2_morethan1,SOMATICSNIPER_morethan1,VARSCAN2_morethan1)
  #Prepare CONSENSUS2 MAF file
  Total.dt <- data.table(Total_morethan1)
  setkey(Total.dt, name)
  Total.dt$t_depth<-as.numeric(as.character(Total.dt$t_depth))
  Total.dt$t_ref_count<-as.numeric(as.character(Total.dt$t_ref_count))
  Total.dt$t_alt_count<-as.numeric(as.character(Total.dt$t_alt_count))
  Total.dt$n_depth<-as.numeric(as.character(Total.dt$n_depth))
  #Define mean reference and variant alleles read depths
  CONSENSUS2<-data.frame(Total.dt[, c(list(Hugo_Symbol=Hugo_Symbol[1],Entrez_Gene_Id=Entrez_Gene_Id[1],Center=Center[1],NCBI_Build=NCBI_Build[1],Chromosome=Chromosome[1],Start_Position=Start_Position[1],End_Position=End_Position[1],Strand=Strand[1],Variant_Classification=Variant_Classification[1],Variant_Type=Variant_Type[1],Reference_Allele=Reference_Allele[1],Tumor_Seq_Allele1=Tumor_Seq_Allele1[1],Tumor_Seq_Allele2=Tumor_Seq_Allele2[1],dbSNP_RS=dbSNP_RS[1],dbSNP_Val_Status=dbSNP_Val_Status[1],Tumor_Sample_Barcode=Tumor_Sample_Barcode[1],Matched_Norm_Sample_Barcode=Matched_Norm_Sample_Barcode[1],Match_Norm_Seq_Allele1=Match_Norm_Seq_Allele1[1],Match_Norm_Seq_Allele2=Match_Norm_Seq_Allele2[1],Tumor_Validation_Allele1=Tumor_Validation_Allele1[1],Tumor_Validation_Allele2=Tumor_Validation_Allele2[1],Match_Norm_Validation_Allele1=Match_Norm_Validation_Allele1[1],Match_Norm_Validation_Allele2=Match_Norm_Validation_Allele2[1],Verification_Status=Verification_Status[1],Validation_Status=Validation_Status[1],Mutation_Status=Mutation_Status[1],Sequencing_Phase=Sequencing_Phase[1],Sequence_Source=Sequence_Source[1],Validation_Method=Validation_Method[1],Score=Score[1],BAM_File=BAM_File[1],Sequencer=Sequencer[1],Tumor_Sample_UUID=Tumor_Sample_UUID[1],Matched_Norm_Sample_UUID=Matched_Norm_Sample_UUID[1],HGVSc=HGVSc[1],HGVSp=HGVSp[1],HGVSp_Short=HGVSp_Short[1],Transcript_ID=Transcript_ID[1],Exon_Number=Exon_Number[1],n_ref_count=n_ref_count[1],n_alt_count=n_alt_count[1],all_effects=all_effects[1],Allele=Allele[1],Gene=Gene[1],Feature=Feature[1],Feature_type=Feature_type[1],One_Consequence=One_Consequence[1],Consequence=Consequence[1],cDNA_position=cDNA_position[1],CDS_position=CDS_position[1],Protein_position=Protein_position[1],Amino_acids=Amino_acids[1],Codons=Codons[1],Existing_variation=Existing_variation[1],ALLELE_NUM=ALLELE_NUM[1],DISTANCE=DISTANCE[1],TRANSCRIPT_STRAND=TRANSCRIPT_STRAND[1],SYMBOL=SYMBOL[1],SYMBOL_SOURCE=SYMBOL_SOURCE[1],HGNC_ID=HGNC_ID[1],BIOTYPE=BIOTYPE[1],CANONICAL=CANONICAL[1],CCDS=CCDS[1],ENSP=ENSP[1],SWISSPROT=SWISSPROT[1],TREMBL=TREMBL[1],UNIPARC=UNIPARC[1],RefSeq=RefSeq[1],SIFT=SIFT[1],PolyPhen=PolyPhen[1],EXON=EXON[1],INTRON=INTRON[1],DOMAINS=DOMAINS[1],GMAF=GMAF[1],AFR_MAF=AFR_MAF[1],AMR_MAF=AMR_MAF[1],ASN_MAF=ASN_MAF[1],EAS_MAF=EAS_MAF[1],EUR_MAF=EUR_MAF[1],SAS_MAF=SAS_MAF[1],AA_MAF=AA_MAF[1],EA_MAF=EA_MAF[1],CLIN_SIG=CLIN_SIG[1],SOMATIC=SOMATIC[1],PUBMED=PUBMED[1],MOTIF_NAME=MOTIF_NAME[1],MOTIF_POS=MOTIF_POS[1],HIGH_INF_POS=HIGH_INF_POS[1],MOTIF_SCORE_CHANGE=MOTIF_SCORE_CHANGE[1],IMPACT=IMPACT[1],PICK=PICK[1],VARIANT_CLASS=VARIANT_CLASS[1],TSL=TSL[1],HGVS_OFFSET=HGVS_OFFSET[1],PHENO=PHENO[1],MINIMISED=MINIMISED[1],ExAC_AF=ExAC_AF[1],ExAC_AF_Adj=ExAC_AF_Adj[1],ExAC_AF_AFR=ExAC_AF_AFR[1],ExAC_AF_AMR=ExAC_AF_AMR[1],ExAC_AF_EAS=ExAC_AF_EAS[1],ExAC_AF_FIN=ExAC_AF_FIN[1],ExAC_AF_NFE=ExAC_AF_NFE[1],ExAC_AF_OTH=ExAC_AF_OTH[1],ExAC_AF_SAS=ExAC_AF_SAS[1],GENE_PHENO=GENE_PHENO[1],FILTER=FILTER[1],CONTEXT=CONTEXT[1],src_vcf_id=src_vcf_id[1],tumor_bam_uuid=tumor_bam_uuid[1],normal_bam_uuid=normal_bam_uuid[1],case_id=case_id[1],GDC_FILTER=GDC_FILTER[1],COSMIC=COSMIC[1],MC3_Overlap=MC3_Overlap[1],GDC_Validation_Status=GDC_Validation_Status[1],sum=sum[1]), lapply(.SD,mean)), by ="name", .SDcols=c("t_depth","t_ref_count","t_alt_count","n_depth")])
  CONSENSUS2$t_depth<-round(CONSENSUS2$t_depth)
  CONSENSUS2$t_ref_count<-round(CONSENSUS2$t_ref_count)
  CONSENSUS2$t_alt_count<-round(CONSENSUS2$t_alt_count)
  CONSENSUS2$n_depth<-round(CONSENSUS2$n_depth)
  #Final CONSENSUS2 MAF file
  CONSENSUS2<-CONSENSUS2[,c(1:40,119:122,41:118)]
  CONSENSUS2<-CONSENSUS2[,-c(1,122)]
  CONSENSUS2<-CONSENSUS2[order(CONSENSUS2$Chromosome, CONSENSUS2$Start_Position), ]
  #Write CONSENSUS2 MAF file
  write.table(CONSENSUS2, paste0(dataDir,"/Maf_files/TCGA_",cancer[i],"_CONSENSUS2.maf"), row.names = F, col.names = T, quote = F,sep="\t")
}

################################################################################
### Step 5: Create UNION MAF file with all variants detected by any variant caller 
################################################################################
print("Create UNION MAF file with all variants detected by any variant caller")
## Read variants detected by all variant callers 
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  for(k in 1:length(caller)){
    print(paste(i,cancer[i],caller[k],sep = "_"))
    maf<-NULL
    maf_union<-NULL
    filename<-NULL
    union<-NULL
    #Set filename
    filename<-paste(cancer[i],caller[k],sep=".")
    #Read maf
    maf<-read.maf(paste0(dataDir,"/","Maf_files","/","TCGA","_",cancer[i],"_",caller[k],".maf.gz"), vc_nonSyn=c("Splice_Region","Intron","5'Flank","Silent","3'UTR","RNA","5'UTR","IGR","3'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site","Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", "Missense_Mutation"),removeDuplicatedVariants=F)
    #Create unique variant ID
    maf@data$name<-paste(maf@data$Chromosome, maf@data$Start_Position, maf@data$End_Position, maf@data$Tumor_Seq_Allele1, maf@data$Tumor_Seq_Allele2, maf@data$Tumor_Sample_Barcode, sep = ":")
    #Format data
    maf@data$caller<-1
    colnames(maf@data)[122]<-paste0(caller[k])
    #Select Union of all variant callers
    union<-read.table(paste0(dataDir,"/Maf_files/",cancer[i],"_comparison.txt"), header = T, sep = "\t")
    maf_union<-merge(maf@data, union[,c(1,6)], by = "name")
    #Format data
    maf_union$src_vcf_id<-"NA"
    maf_union$Caller<-paste0(caller[k])
    maf_union<-maf_union[,-122]
    #Write results
    write.table(maf_union, paste0(dataDir,"/Maf_files/",filename,".Union.txt"), row.names = F, col.names = T, quote = F,sep="\t")
  }
}

## Prepare UNION MAF file
for(i in 1:length(cancer)){
  #print(i)
  #print(cancer[i])
  MUSE_union<-NULL
  MUTECT2_union<-NULL
  SOMATICSNIPER_union<-NULL
  VARSCAN2_union<-NULL
  Total_union<-NULL
  Total.dt<-NULL
  UNION<-NULL
  #Read variants detected by each caller that are also detected by at least two other variant callers
  MUSE_union<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUSE.Union.txt"))
  MUTECT2_union<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".MUTECT2.Union.txt"))
  SOMATICSNIPER_union<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".SOMATICSNIPER.Union.txt"))
  VARSCAN2_union<-data.table::fread(paste0(dataDir,"/Maf_files/",cancer[i],".VARSCAN2.Union.txt"))
  #Combine all variants detected by all callers
  Total_union<-rbind(MUSE_union,MUTECT2_union,SOMATICSNIPER_union,VARSCAN2_union)
  #Prepare UNION MAF file
  Total.dt <- data.table(Total_union)
  setkey(Total.dt, name)
  Total.dt$t_depth<-as.numeric(as.character(Total.dt$t_depth))
  Total.dt$t_ref_count<-as.numeric(as.character(Total.dt$t_ref_count))
  Total.dt$t_alt_count<-as.numeric(as.character(Total.dt$t_alt_count))
  Total.dt$n_depth<-as.numeric(as.character(Total.dt$n_depth))
  #Define mean reference and variant alleles read depths
  UNION<-data.frame(Total.dt[, c(list(Hugo_Symbol=Hugo_Symbol[1],Entrez_Gene_Id=Entrez_Gene_Id[1],Center=Center[1],NCBI_Build=NCBI_Build[1],Chromosome=Chromosome[1],Start_Position=Start_Position[1],End_Position=End_Position[1],Strand=Strand[1],Variant_Classification=Variant_Classification[1],Variant_Type=Variant_Type[1],Reference_Allele=Reference_Allele[1],Tumor_Seq_Allele1=Tumor_Seq_Allele1[1],Tumor_Seq_Allele2=Tumor_Seq_Allele2[1],dbSNP_RS=dbSNP_RS[1],dbSNP_Val_Status=dbSNP_Val_Status[1],Tumor_Sample_Barcode=Tumor_Sample_Barcode[1],Matched_Norm_Sample_Barcode=Matched_Norm_Sample_Barcode[1],Match_Norm_Seq_Allele1=Match_Norm_Seq_Allele1[1],Match_Norm_Seq_Allele2=Match_Norm_Seq_Allele2[1],Tumor_Validation_Allele1=Tumor_Validation_Allele1[1],Tumor_Validation_Allele2=Tumor_Validation_Allele2[1],Match_Norm_Validation_Allele1=Match_Norm_Validation_Allele1[1],Match_Norm_Validation_Allele2=Match_Norm_Validation_Allele2[1],Verification_Status=Verification_Status[1],Validation_Status=Validation_Status[1],Mutation_Status=Mutation_Status[1],Sequencing_Phase=Sequencing_Phase[1],Sequence_Source=Sequence_Source[1],Validation_Method=Validation_Method[1],Score=Score[1],BAM_File=BAM_File[1],Sequencer=Sequencer[1],Tumor_Sample_UUID=Tumor_Sample_UUID[1],Matched_Norm_Sample_UUID=Matched_Norm_Sample_UUID[1],HGVSc=HGVSc[1],HGVSp=HGVSp[1],HGVSp_Short=HGVSp_Short[1],Transcript_ID=Transcript_ID[1],Exon_Number=Exon_Number[1],n_ref_count=n_ref_count[1],n_alt_count=n_alt_count[1],all_effects=all_effects[1],Allele=Allele[1],Gene=Gene[1],Feature=Feature[1],Feature_type=Feature_type[1],One_Consequence=One_Consequence[1],Consequence=Consequence[1],cDNA_position=cDNA_position[1],CDS_position=CDS_position[1],Protein_position=Protein_position[1],Amino_acids=Amino_acids[1],Codons=Codons[1],Existing_variation=Existing_variation[1],ALLELE_NUM=ALLELE_NUM[1],DISTANCE=DISTANCE[1],TRANSCRIPT_STRAND=TRANSCRIPT_STRAND[1],SYMBOL=SYMBOL[1],SYMBOL_SOURCE=SYMBOL_SOURCE[1],HGNC_ID=HGNC_ID[1],BIOTYPE=BIOTYPE[1],CANONICAL=CANONICAL[1],CCDS=CCDS[1],ENSP=ENSP[1],SWISSPROT=SWISSPROT[1],TREMBL=TREMBL[1],UNIPARC=UNIPARC[1],RefSeq=RefSeq[1],SIFT=SIFT[1],PolyPhen=PolyPhen[1],EXON=EXON[1],INTRON=INTRON[1],DOMAINS=DOMAINS[1],GMAF=GMAF[1],AFR_MAF=AFR_MAF[1],AMR_MAF=AMR_MAF[1],ASN_MAF=ASN_MAF[1],EAS_MAF=EAS_MAF[1],EUR_MAF=EUR_MAF[1],SAS_MAF=SAS_MAF[1],AA_MAF=AA_MAF[1],EA_MAF=EA_MAF[1],CLIN_SIG=CLIN_SIG[1],SOMATIC=SOMATIC[1],PUBMED=PUBMED[1],MOTIF_NAME=MOTIF_NAME[1],MOTIF_POS=MOTIF_POS[1],HIGH_INF_POS=HIGH_INF_POS[1],MOTIF_SCORE_CHANGE=MOTIF_SCORE_CHANGE[1],IMPACT=IMPACT[1],PICK=PICK[1],VARIANT_CLASS=VARIANT_CLASS[1],TSL=TSL[1],HGVS_OFFSET=HGVS_OFFSET[1],PHENO=PHENO[1],MINIMISED=MINIMISED[1],ExAC_AF=ExAC_AF[1],ExAC_AF_Adj=ExAC_AF_Adj[1],ExAC_AF_AFR=ExAC_AF_AFR[1],ExAC_AF_AMR=ExAC_AF_AMR[1],ExAC_AF_EAS=ExAC_AF_EAS[1],ExAC_AF_FIN=ExAC_AF_FIN[1],ExAC_AF_NFE=ExAC_AF_NFE[1],ExAC_AF_OTH=ExAC_AF_OTH[1],ExAC_AF_SAS=ExAC_AF_SAS[1],GENE_PHENO=GENE_PHENO[1],FILTER=FILTER[1],CONTEXT=CONTEXT[1],src_vcf_id=src_vcf_id[1],tumor_bam_uuid=tumor_bam_uuid[1],normal_bam_uuid=normal_bam_uuid[1],case_id=case_id[1],GDC_FILTER=GDC_FILTER[1],COSMIC=COSMIC[1],MC3_Overlap=MC3_Overlap[1],GDC_Validation_Status=GDC_Validation_Status[1],sum=sum[1]), lapply(.SD,mean)), by ="name", .SDcols=c("t_depth","t_ref_count","t_alt_count","n_depth")])
  UNION$t_depth<-round(UNION$t_depth)
  UNION$t_ref_count<-round(UNION$t_ref_count)
  UNION$t_alt_count<-round(UNION$t_alt_count)
  UNION$n_depth<-round(UNION$n_depth)
  #Final UNION MAF file
  UNION<-UNION[,c(1:40,119:122,41:118)]
  UNION<-UNION[,-c(1,122)]
  UNION<-UNION[order(UNION$Chromosome, UNION$Start_Position), ]
  #Write UNION MAF file
  write.table(UNION, paste0(dataDir,"/Maf_files/TCGA_",cancer[i],"_UNION.maf"), row.names = F, col.names = T, quote = F,sep="\t")
}

################################################################################
### Step 6: Remove intermediate files 
################################################################################
print("Removing intermediate files")
# Remove intermediate files
unlink(paste0(dataDir,"/Maf_files/*txt"))
print("¡¡¡Analysis finished!!!")
#===============================================================================