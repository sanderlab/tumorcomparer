# Library
library(tidyverse)

# Working directory 
setwd("~/Google Drive/tumorcomparer/pancan_data/tumorcomparer_prepare_scripts/")
stopifnot(file.exists("cclp_annotated_maf.txt"))

# Load data
cclp_cosmic_annovar <- read_tsv("cclp_annotated_maf.txt")
tcga_cosmic_annovar <- read_tsv("tcga_annotated_maf.txt")

# CCLP
cclp_cosmic_annovar_sm <- cclp_cosmic_annovar[which(cclp_cosmic_annovar$TumorComparer_Binary_Status), c("Gene name", "Sample name", "TumorComparer_Binary_Status")]
cclp_genes <- unique(cclp_cosmic_annovar_sm$`Gene name`)
cclp_samples <- unique(cclp_cosmic_annovar_sm$`Sample name`)

cclp_bin_mat <- matrix(0, nrow=length(cclp_genes), ncol=length(cclp_samples), dimnames=list(cclp_genes, cclp_samples))

max_iter <- nrow(cclp_cosmic_annovar_sm)
pb <- txtProgressBar(min=1, max=max_iter, style=3)

for(i in 1:nrow(cclp_cosmic_annovar_sm)) {
  setTxtProgressBar(pb, i)
  
  gene <- as.character(cclp_cosmic_annovar_sm[i, "Gene name"])
  sample <- as.character(cclp_cosmic_annovar_sm[i, "Sample name"])
  
  gene_idx <- which(gene == cclp_genes)
  sample_idx <- which(sample == cclp_samples)
  
  cclp_bin_mat[gene_idx, sample_idx] <- 1
}

cclp_bin_df <- data.frame(cclp_bin_mat, check.names=FALSE)

write.table(cclp_bin_df, "cclp_binary_mutation_matrix.txt", sep="\t", quote=FALSE, row.names=TRUE,
            col.names=TRUE)

## Checks
cclp_bin_df[1:10, 1:10]
a <- length(which(cclp_bin_df == 1))
b <- length(which(cclp_bin_df == 0))

# TCGA 
tcga_cosmic_annovar_sm <- tcga_cosmic_annovar[which(tcga_cosmic_annovar$TumorComparer_Binary_Status), c("Hugo_Symbol", "Tumor_Sample_Barcode", "TumorComparer_Binary_Status")]
tcga_genes <- unique(tcga_cosmic_annovar_sm$`Hugo_Symbol`)
tcga_samples <- unique(tcga_cosmic_annovar_sm$`Tumor_Sample_Barcode`)

tcga_bin_mat <- matrix(0, nrow=length(tcga_genes), ncol=length(tcga_samples), dimnames=list(tcga_genes, tcga_samples))

max_iter <- nrow(tcga_cosmic_annovar_sm)
pb <- txtProgressBar(min=1, max=max_iter, style=3)

for(i in 1:nrow(tcga_cosmic_annovar_sm)) {
  #cat("I: ", i, "\n")
  setTxtProgressBar(pb, i)
  
  gene <- as.character(tcga_cosmic_annovar_sm[i, "Hugo_Symbol"])
  sample <- as.character(tcga_cosmic_annovar_sm[i, "Tumor_Sample_Barcode"])
  
  gene_idx <- which(gene == tcga_genes)
  sample_idx <- which(sample == tcga_samples)
  
  tcga_bin_mat[gene_idx, sample_idx] <- 1
}

tcga_bin_df <- data.frame(tcga_bin_mat, check.names=FALSE)

write.table(tcga_bin_df, "tcga_binary_mutation_matrix.txt", sep="\t", quote=FALSE, row.names=TRUE,
            col.names=TRUE)

## Checks
tcga_bin_df[1:10, 1:10]
a <- length(which(tcga_bin_df == 1))
b <- length(which(tcga_bin_df == 0))

