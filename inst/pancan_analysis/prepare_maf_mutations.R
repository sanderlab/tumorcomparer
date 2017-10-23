# Load libraries
library(tidyverse)
#options(tibble.width = Inf)

# Set working directory
setwd("~/Google Drive/tumorcomparer/pancan_data")

# Parameters 
tcga_mutations_remove <- c("3'Flank", "3'UTR", "5'Flank", "5'UTR", "Intron", "RNA", "Silent")
tcga_mutations_deleterious <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
                                "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")
tcga_mutations_missense <- c("Missense_Mutation")

cclp_mutations_remove <- c("Complex - compound substitution", "Substitution - coding silent", "Unknown")
cclp_mutations_deleterious <- c("Complex - deletion inframe", "Complex - frameshift", "Complex - insertion inframe", 
                                "Deletion - Frameshift", "Deletion - In frame", "Insertion - Frameshift", "Insertion - In frame",
                                "Nonstop extension", "Substitution - Nonsense")
cclp_mutations_missense <- c("Substitution - Missense")

# Load data
filtered_cosmic_mutations <- readRDS("tumorcomparer_prepare_scripts/filtered_cosmic_mutations.rds")

if(!file.exists("tumorcomparer_prepare_scripts/tcga_maf.rds")) {
  tcga <- read_tsv("mc3.v0.2.8.PUBLIC.maf", col_names = TRUE, na = ".")
  
  idx <- which(!(tcga$Variant_Classification %in% tcga_mutations_remove))
  tcga <- tcga[idx, ]
  
  saveRDS(tcga, "tumorcomparer_prepare_scripts/tcga_maf.rds")
} else {
  tcga <- readRDS("tumorcomparer_prepare_scripts/tcga_maf.rds")
}

if(!file.exists("tumorcomparer_prepare_scripts/cclp_maf.rds")) {
  cclp <- read_tsv("CCLP_v81/CosmicCLP_MutantExport.tsv", col_names = TRUE, na = "")
  
  idx <- which(!(cclp$`Mutation Description` %in% cclp_mutations_remove))
  cclp <- cclp[idx, ]

  saveRDS(cclp, "tumorcomparer_prepare_scripts/cclp_maf.rds")
} else {
  cclp <- readRDS("tumorcomparer_prepare_scripts/cclp_maf.rds")
}

if(!file.exists("tumorcomparer_prepare_scripts/cclp_annovar.rds")) {
  cclp_annovar <- read_tsv("annovar_results/gdsc_msmut_annovar.hg19_multianno.txt", 
                          col_names = TRUE, col_types=cols(Chr = "c"), na = ".")
  saveRDS(cclp_annovar, "tumorcomparer_prepare_scripts/cclp_annovar.rds")
} else {
  cclp_annovar <- readRDS("tumorcomparer_prepare_scripts/cclp_annovar.rds")
}

if(!file.exists("tumorcomparer_prepare_scripts/tcga_annovar.rds")) {
  tcga_annovar <- read_tsv("annovar_results/tcga_msmut_annovar.hg19_multianno.txt", 
                          col_names = TRUE, col_types=cols(Chr = "c"), na = ".")
  saveRDS(tcga_annovar, "tumorcomparer_prepare_scripts/tcga_annovar.rds")
} else {
  tcga_annovar <- readRDS("tumorcomparer_prepare_scripts/tcga_annovar.rds")
}

# Filter data 
#tcga <- tcga[c(tcga_mutations_deleterious)]
#cclp <- cclp[c(cclp_mutations_deleterious)]

# Unique samples/genes 
tmp_samples <- unique(tcga$Tumor_Sample_Barcode)
length(tmp_samples)

tmp_genes <- unique(tcga$Hugo_Symbol)
length(tmp_genes)

original_tcga <- tcga

# # Iterate over MAF contents  
# i <- 1
# tmp_sample <- tmp_samples[i]
# sample_idx <- which(tcga$Tumor_Sample_Barcode == tmp_sample)
# tmp <- tcga[sample_idx,]
# tmp_genes <- tmp$Hugo_Symbol
# 
# j <- 1 
# tmp_gene <- tmp_genes[j]
# gene_idx <- which(tcga_annovar$Gene.refGene == tmp_gene)
# tcga_annovar[gene_idx, ]
# 
# tmp_annovar <- tcga_annovar[gene_idx, ]
# ma_high_idx <- which(tmp_annovar$MutationAssessor_pred == "H")
# tmp_ma_high <- tmp_annovar[ma_high_idx, ]
# 
# # Number of Mutation Assessor "High" entries for gene for patient 
# nrow(tmp_ma_high) 

# Process TCGA
#tcga_mutations_missense_idx <- which(tcga$Variant_Classification %in% tcga_mutations_missense)
#tcga_mutations_missense_full <- tcga[tcga_mutations_missense_idx,]

#i <- 1
#tcga_mutations_missense_full[i, c("Hugo_Symbol", "HGVSp_Short")]

# Process CCLP
system.time({
  cclp_cosmic <- merge(cclp, filtered_cosmic_mutations, 
              by.x = c("Gene name", "Mutation AA"), 
              by.y = c("Gene name", "Mutation AA"),
              all.x = TRUE)
})

cclp_minimal_annovar <- cclp_annovar[, c("Chr", "Start", "End", "Gene.refGene", 
                                         "MutationAssessor_score",
                                         "MutationAssessor_pred",
                                         "AAChange.refGene")]

## Grab the amino acid change
tmp_aa <- strsplit(cclp_minimal_annovar$AAChange.refGene, ":")
tmp_aa2 <- lapply(tmp_aa, function(x) {
  x[length(x)]
})
tmp_aa3 <- unlist(tmp_aa2)
cclp_minimal_annovar$AAChange <- tmp_aa3

system.time({
  cclp_cosmic_annovar <- merge(cclp_cosmic, cclp_minimal_annovar, 
              by.x = c( "Gene name", "Mutation AA"), 
              by.y = c("Gene.refGene", "AAChange"),
              all.x = TRUE)
})

cclp_cosmic_annovar$TumorComparer_Binary_Status <- with(cclp_cosmic_annovar, ifelse(`Mutation Description` %in% cclp_mutations_deleterious |
                                                    (!is.na(COSMIC_StudyCount) & COSMIC_StudyCount >= 10) |
                                                    (!is.na(MutationAssessor_pred) & MutationAssessor_pred == "H"),
                                                  TRUE, FALSE))
# DEBUG
# m2 <- cclp_cosmic_annovar[which(cclp_cosmic_annovar$TumorComparer_Binary_Status), c("Gene name", "Mutation AA", "Mutation Description", "COSMIC_StudyCount", 
#                                                   "MutationAssessor_pred", "TumorComparer_Binary_Status")]
# 
# m3 <- cclp_cosmic_annovar[, c("Gene name", "Mutation AA", "Mutation Description", "COSMIC_StudyCount", 
#              "MutationAssessor_pred", "TumorComparer_Binary_Status")]
# 
# m4 <- cclp_cosmic_annovar[with(z, order(`Mutation Description`)), c("Gene name", "Mutation AA", "Mutation Description", "COSMIC_StudyCount", 
#                                                  "MutationAssessor_pred", "TumorComparer_Binary_Status")]

write_tsv(cclp_cosmic_annovar, "tumorcomparer_prepare_scripts/cclp_annotated_maf.txt")

# Process TCGA
system.time({
  tcga_cosmic <- merge(tcga, filtered_cosmic_mutations, 
        by.x = c("Hugo_Symbol", "HGVSp_Short"), 
        by.y = c("Gene name", "Mutation AA"),
        all.x = TRUE)
})

minimal_annovar <- tcga_annovar[, c("Chr", "Start", "End", "Gene.refGene", 
                                    "MutationAssessor_score",
                                    "MutationAssessor_pred")]

system.time({
  tcga_cosmic_annovar <- merge(tcga_cosmic, minimal_annovar, 
             by.x = c("Chromosome", "Start_Position", "End_Position", "Hugo_Symbol"), 
             by.y = c("Chr", "Start", "End", "Gene.refGene"),
             all.x = TRUE)
})

tcga_cosmic_annovar$TumorComparer_Binary_Status <- with(tcga_cosmic_annovar, ifelse(Variant_Classification %in% tcga_mutations_deleterious |
                                                (!is.na(COSMIC_StudyCount) & COSMIC_StudyCount >= 10) |
                                                (!is.na(MutationAssessor_pred) & MutationAssessor_pred == "H"),
                                                TRUE, FALSE))

# DEBUG
# z2 <- z[which(z$TumorComparer_Binary_Status), c("Hugo_Symbol", "HGVSp_Short", "Variant_Classification", "COSMIC_StudyCount", 
#                                                 "MutationAssessor_pred", "TumorComparer_Binary_Status")]
# 
# z3 <- z[, c("Hugo_Symbol", "HGVSp_Short", "Variant_Classification", "COSMIC_StudyCount", 
#             "MutationAssessor_pred", "TumorComparer_Binary_Status")]
# 
# z4 <- z[with(z, order(Variant_Classification)), c("Hugo_Symbol", "HGVSp_Short", "Variant_Classification", "COSMIC_StudyCount", 
#                                                   "MutationAssessor_pred", "TumorComparer_Binary_Status")]

write_tsv(tcga_cosmic_annovar, "tumorcomparer_prepare_scripts/tcga_annotated_maf.txt")

