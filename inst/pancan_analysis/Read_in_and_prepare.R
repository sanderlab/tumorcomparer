#library(stringr)
library(matrixStats)
#library(preprocessCore)

# SET WORKING DIRECTORY ----
setwd("~/Google\ Drive/tumorcomparer/pancan_data/")

# LOAD MISC FUNCTIONS ----
system.file("pancan_analysis", "tumorcomparer_prepare_scripts/misc_functions.R", package = "tumorcomparer")

# READ IN DATA ----
TCGA_PanCan_RNASeq_RPKMs_Weights <- read.delim("tumorcomparer_prepare_scripts/TCGA_PanCan_RNASeq_RPKMs_Weights.txt",row.names = 1)
mut_mat_TCGA_after_Annovar <- read.delim("tumorcomparer_prepare_scripts/tcga_binary_mutation_matrix.txt",check.names = F,row.names = 1)
colnames(mut_mat_TCGA_after_Annovar) <- sapply(colnames(mut_mat_TCGA_after_Annovar),return_3_field_TCGA_ID)
gc()
TCGA_Pancan_all_thresholded_by_genes_whitelisted <- read.delim("all_thresholded.by_genes_whitelisted.tsv",header = T,row.names = 1,check.names = F)
gc()
TCGA_PanCan_annotations <- read.delim("merged_sample_quality_annotations.tsv",check.names = F,header=T)
gc()
TCGA_PanCan_RNASeq_RPKMs <- read.delim("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",header = T,row.names = 1,check.names = F)
gc()
gc()
CCLP_Expression_Microarray <- read.csv("CCLP_v81/CCLP_Gene_Expression_1018_Cell_Lines_19549_Genes.csv",header = T,check.names = F)
CCLP_Expression_Microarray <- CCLP_Expression_Microarray[-14820,]
rownames(CCLP_Expression_Microarray) <- CCLP_Expression_Microarray$Gene
CCLP_Expression_Microarray <- CCLP_Expression_Microarray[,-1]
gc()
CCLP_Cell_Lines_And_Tissue_Types <- read.delim("CCLP_Cell_Lines_And_Tissue_Types.txt",check.names = F)
gc()
mut_mat_CCLP_after_Annovar <- read.delim("tumorcomparer_prepare_scripts/cclp_binary_mutation_matrix.txt",check.names = F,row.names = 1)
gc()
CCLP_GISTIC_thresholded <- read.delim(file="GISTIC/cclp/CCLP_log2_Total_CN_divided_by_2_no_XorY.all_thresholded.by_genes.txt",check.names = F,row.names = 1,header = T)
gc()

# EXTRACT TCGA IDS
TCGA_id_and_tumor_type <- cbind(as.character(TCGA_PanCan_annotations$patient_barcode),as.character(TCGA_PanCan_annotations$`cancer type`))[which(TCGA_PanCan_annotations$Do_not_use == "False"),]
TCGA_id_and_tumor_type <- TCGA_id_and_tumor_type[which(!duplicated(TCGA_id_and_tumor_type[,1])),]
rownames(TCGA_id_and_tumor_type) <- TCGA_id_and_tumor_type[,1]
TCGA_normal_IDs <- colnames(TCGA_PanCan_RNASeq_RPKMs)[grep("11",sapply(colnames(TCGA_PanCan_RNASeq_RPKMs),return_fourth_part_of_TCGA_ID))]
TCGA_tumor_IDs <- setdiff(colnames(TCGA_PanCan_RNASeq_RPKMs),TCGA_normal_IDs)


rownames(TCGA_PanCan_RNASeq_RPKMs)[16302] <- "SLC35E2B|728661"
TCGA_PanCan_RNASeq_RPKMs <- TCGA_PanCan_RNASeq_RPKMs[-c(1:29),]
rownames(TCGA_PanCan_RNASeq_RPKMs) <- sapply(rownames(TCGA_PanCan_RNASeq_RPKMs),return_first_part_Bar)
NAs_per_gene_TCGA <- apply(TCGA_PanCan_RNASeq_RPKMs,1,function(x){length(which(is.na(x)))})
TCGA_PanCan_RNASeq_RPKMs <- TCGA_PanCan_RNASeq_RPKMs[which(NAs_per_gene_TCGA == 0),]
#TCGA_PanCan_RNASeq_RPKMs_Tumors_Zscores <- t(scale(t(TCGA_PanCan_RNASeq_RPKMs[,intersect(colnames(TCGA_PanCan_RNASeq_RPKMs),TCGA_tumor_IDs)])))
TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Zscores <- t(scale(t(log2(1+TCGA_PanCan_RNASeq_RPKMs[,intersect(colnames(TCGA_PanCan_RNASeq_RPKMs),TCGA_tumor_IDs)]))))

# PROCESS CCLP EXPRESSION DATA ----
## z-score data and discretize
system.time(CCLP_Expression_Microarray_Zscores <- t(scale(t(CCLP_Expression_Microarray))))
CCLP_Expression_Microarray_Discretized_By_Zscores <- CCLP_Expression_Microarray_Zscores
CCLP_Expression_Microarray_Discretized_By_Zscores[which(abs(CCLP_Expression_Microarray_Discretized_By_Zscores) < 2)] <- 0
CCLP_Expression_Microarray_Discretized_By_Zscores[which(CCLP_Expression_Microarray_Discretized_By_Zscores < -2)] <- -1
CCLP_Expression_Microarray_Discretized_By_Zscores[which(CCLP_Expression_Microarray_Discretized_By_Zscores > 2)] <- 1

# PROCESS TCGA EXPRESSION DATA ----
TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores <- TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Zscores
TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores[which(abs(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores) < 2)] <- 0
TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores[which(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores < -2)] <- -1
TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores[which(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores > 2)] <- 1
colnames(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores) <- sapply(colnames(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores),return_3_field_TCGA_ID)

# QUANTILE NORMALIZE DATA ----
CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized <- normalize.quantiles(cbind(as.matrix(CCLP_Expression_Microarray[intersect(rownames(CCLP_Expression_Microarray),rownames(TCGA_PanCan_RNASeq_RPKMs[,TCGA_tumor_IDs])),]), as.matrix(log2(1+TCGA_PanCan_RNASeq_RPKMs[intersect(rownames(CCLP_Expression_Microarray),rownames(TCGA_PanCan_RNASeq_RPKMs)),TCGA_tumor_IDs]))))
rownames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized) <- intersect(rownames(CCLP_Expression_Microarray),rownames(TCGA_PanCan_RNASeq_RPKMs[,TCGA_tumor_IDs]))
colnames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized) <- c(colnames(CCLP_Expression_Microarray),sapply(TCGA_tumor_IDs,return_3_field_TCGA_ID))
length(colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted)[grep("11",sapply(colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),return_fourth_part_of_TCGA_ID))]) # confirm that there are  no normals among the CNA samples
colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted) <- sapply(colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),return_3_field_TCGA_ID)

# GET ASSIGNMENT OF CELL LINE/TUMOR ----
GDSC_1000_Binary_Event_Matrix <- read.csv("tumorcomparer_prepare_scripts/GDSC_1000_Binary_Event_Matrix.csv",check.names = F,row.names = 1)
## Get cancer list for all represented cancer types
cancers_list <- intersect(names(table(as.matrix(GDSC_1000_Binary_Event_Matrix)[1,which(GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])),names(table(TCGA_id_and_tumor_type[,2])))
cancers_list <- sort(c(cancers_list,"COAD","READ")) # Because CCLP has COAD/READ together, while TCGA PanCan treats them separately
cancers_list <- cancers_list[-which(cancers_list == "ACC")] #only one ACC cell line in CCLP

# GET GENE SETS ----
COSMIC_CGC_Feb_2017 <- read.csv("tumorcomparer_prepare_scripts/Census_allTue_Feb_21_2017.csv")
Human_Kinase_List <- read.csv("tumorcomparer_prepare_scripts/Human_Kinase_List.csv")
Vogelstein_Drivers <- read.delim("tumorcomparer_prepare_scripts/Driver_Genes_Affected_By_Subtle_Mutations.tsv")

gene_set_exp <- rownames(TCGA_PanCan_RNASeq_RPKMs)[which(rowVars(as.matrix(log2(1+TCGA_PanCan_RNASeq_RPKMs))) > 5)]
tcga_mut_freq <- apply(mut_mat_TCGA_after_Annovar,1,function(x){length(which(x != 0))/length(x)})
cclp_mut_freq <- apply(mut_mat_CCLP_after_Annovar,1,function(x){length(which(x != 0))/length(x)})

freq_exp_dysreg_CCLP <- apply(CCLP_Expression_Microarray_Discretized_By_Zscores,1,function(x){length(which(x!=0))/length(x)})
freq_exp_dysreg_TCGA <- apply(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores,1,function(x){length(which(x!=0))/length(x)})
tcga_cna_amp_freq <- apply(TCGA_Pancan_all_thresholded_by_genes_whitelisted[,-c(1:2)],1,function(x){length(which(x == 2))/length(x)})
tcga_cna_deepdel_freq <- apply(TCGA_Pancan_all_thresholded_by_genes_whitelisted[,-c(1:2)],1,function(x){length(which(x == -2))/length(x)})

gene_set_mut <- rownames(mut_mat_TCGA_after_Annovar)[which(tcga_mut_freq >= 0.01)]
gene_set_cna <- rownames(CCLP_GISTIC_thresholded)[which(tcga_cna_amp_freq >= 0.05 | tcga_cna_deepdel_freq >= 0.02)]
