source("Weighted_similarity_correlation.r")

# The Jaccard measure
#source("Weighted_similarity_own_measure.r")

# Files 
CCLE_..._CNA.txt: For cell lines  

Replace this 

tumor_MUT <- read.table("tumor_MUT.txt",sep="\t",header=T,row.names=1)
tumor_CNA <- read.table("tumor_CNA.txt",sep="\t",header=T,row.names=1)
cell_line_MUT <- read.table("cell_line_MUT.txt",sep="\t",header=T,row.names=1)
cell_line_CNA <- read.table("cell_line_CNA.txt",sep="\t",header=T,row.names=1)

Pan-cancer Weights: Default_weights_for_known_cancer_genes.txt
Specific disease: Genes_and_weights_TCGA_BRCA_based

# Change all these to a specific cancer type 

tumor_MUT.txt
tumor_CNA.txt
cell_line_MUT.txt
cell_line_CNA.txt
Genes_and_weights.txt
Default_weights_for_known_cancer_genes.txt (ALWAYS THE SAME)

# Default weight if no MutSig and GISTIC 
pseudoweight <- 0.000001

# Included in the weighted similarity correlation
K_nearest_neighbors_etc.r

# Dumping the composite_alteration_matrix
composite_alteration_matrix.txt 

# Uses subtypes (but ignored for now)
CCLE_vs_TCGA_Breast_Cancer_Unweighted.png


