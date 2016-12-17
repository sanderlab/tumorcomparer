# LOAD DATA ---- 
CNA_default_weight <- 0.01
MUT_default_weight <- 0.01

CNA_known_cancer_gene_weight <- 0.1
MUT_known_cancer_gene_weight <- 0.1

tumor_mut_file <- "tumor_MUT.txt"
tumor_cna_file <- "tumor_CNA.txt"
cell_line_mut_file <- "cell_line_MUT.txt"
cell_line_cna_file <- "cell_line_CNA.txt"

pancancer_gene_weights_file <- "Default_weights_for_known_cancer_genes.txt"
cancer_specific_gene_weights_file <- "Genes_and_weights.txt"

output_composite_alteration_matrix_file <- "composite_alteration_matrix.txt"

fraction_of_tumors_for_comparison <- 0.1 

# GET INTERSECTING GENE ----
tumor_MUT <- read.table(tumor_mut_file,sep="\t",header=TRUE,row.names=1)
tumor_CNA <- read.table(tumor_cna_file,sep="\t",header=TRUE,row.names=1)
cell_line_MUT <- read.table(cell_line_mut_file,sep="\t",header=TRUE,row.names=1)
cell_line_CNA <- read.table(cell_line_cna_file,sep="\t",header=TRUE,row.names=1)

tumors_with_both_MUT_and_CNA <- intersect(colnames(tumor_MUT),colnames(tumor_CNA))
cell_lines_with_both_MUT_and_CNA <- intersect(colnames(cell_line_MUT),colnames(cell_line_CNA))
cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, return_first_part)

genes_with_MUT_in_both <- intersect(rownames(tumor_MUT),rownames(cell_line_MUT))
genes_with_CNA_in_both <- intersect(rownames(tumor_CNA),rownames(cell_line_CNA))

# Need not do this unless really want MUT and CNA data for the same gene sets
genes_in_all_4_files <- intersect(genes_with_MUT_in_both,genes_with_CNA_in_both) 

cell_line_CNA_high_level_only <- apply(cell_line_CNA,2,keep_only_high_level_cnas)
tumor_CNA_high_level_only <- apply(tumor_CNA,2,keep_only_high_level_cnas)

# CREATE COMPOSITE MATRIX ----
## Rows genes (MUT, CNA), columns (Tumors/cell lines)
composite_CNA <- cbind(cell_line_CNA[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
composite_CNA_high_level_only <- cbind(cell_line_CNA_high_level_only[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA_high_level_only[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
composite_MUT <- cbind(cell_line_MUT[genes_with_MUT_in_both,cell_lines_with_both_MUT_and_CNA],tumor_MUT[genes_with_MUT_in_both,tumors_with_both_MUT_and_CNA])

rownames(composite_MUT) <- paste(rownames(composite_MUT),"MUT",sep="_")
rownames(composite_CNA) <- paste(rownames(composite_CNA),"CNA",sep="_")
rownames(composite_CNA_high_level_only) <- paste(rownames(composite_CNA_high_level_only),"CNA",sep="_")

composite_mat <- rbind(composite_MUT,composite_CNA)
alt_mat <- rbind(composite_MUT,composite_CNA)

# WRITE COMPOSITE
write.table(composite_mat,file=output_composite_alteration_matrix_file,sep="\t",quote=FALSE)
#composite_mat_high_level_only <- rbind(composite_MUT,composite_CNA_high_level_only)

# Calculation of alteration frequencies 
# Assign frequency weights as (freq. of alteration of gene)/(mean freq. of alteration across all genes) - "rewarding recurrent changes"
overall_alt_freq <- length(which((alt_mat[]) != 0))/ ( length(which((alt_mat[])==0)) + length(which((alt_mat[])!=0)))

freq_alt <- rep(0, nrow(alt_mat))
freq_alt <- apply(alt_mat,1,compute_freq_alt)
#freq_alt_high_level  <- apply(composite_mat_high_level_only,1,compute_freq_alt)
freq_alt_mut_tumors <- apply(composite_MUT[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)
freq_alt_cna_tumors <- apply(composite_CNA[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)

freq_alt_samplewise <- apply(alt_mat,2,compute_freq_alt)
#freq_alt_samplewise_CNA_high_level_only <- apply(composite_mat_high_level_only,2,compute_freq_alt)
alt_mat <- alt_mat[,which(freq_alt_samplewise > 0)]
composite_mat <- composite_mat[,which(freq_alt_samplewise > 0)]
#composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise > 0)]

names(freq_alt) <- rownames(alt_mat)
freq_weights <- rep(1,nrow(alt_mat))
freq_weights <- freq_alt/overall_alt_freq
names(freq_weights) <- rownames(alt_mat)


# GET WEIGHTS 
# Read in user-provided weights
known_cancer_genes_and_weights_all <- read.table(pancancer_gene_weights_file, sep="\t",header=TRUE,row.names=1) 
known_cancer_genes_and_weights <- as.matrix(known_cancer_genes_and_weights_all[intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
genes_and_weights_all <- read.table(cancer_specific_gene_weights_file,sep="\t",header=TRUE,row.names=1) 
genes_and_weights <- as.matrix(genes_and_weights_all[intersect(rownames(genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
rownames(genes_and_weights) <- intersect(rownames(genes_and_weights_all),rownames(alt_mat))
rownames(known_cancer_genes_and_weights) <- intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat))

annotation_weights <- rep(NA,nrow(alt_mat)) 
annotation_weights[grep("_CNA",rownames(alt_mat))] <- CNA_default_weight
annotation_weights[grep("_MUT",rownames(alt_mat))] <- MUT_default_weight

names(annotation_weights) <- rownames(alt_mat) 
for(i in 1:nrow(known_cancer_genes_and_weights)) # Overwrite default weights for known cancer genes
  annotation_weights[rownames(known_cancer_genes_and_weights)[i]] = known_cancer_genes_and_weights[i,]
for(i in 1:nrow(genes_and_weights)) # Overwrite default weights for input provided weights
  annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i,]


gene_weights <- rep(1,nrow(alt_mat))
names(gene_weights) <- rownames(alt_mat)

gene_weights <- annotation_weights # if using user-provided weights only

gene_weights <- gene_weights/max(gene_weights) # map to 0-1

# CALCULATE CORRELATIONS ----
# Including low-level CNAs
cor_weighted <- calc_weighted_corr(as.matrix(composite_mat),as.matrix(composite_mat),gene_weights)
# Excluding low-levels CNAs
#cor_weighted_high_level_only <- calc_weighted_corr(as.matrix(composite_mat_high_level_only),as.matrix(composite_mat_high_level_only),gene_weights)
cor_unweighted <- cor(alt_mat)

# CATEGORIZE ----
## Set K-nearest neighbors 
k <- fraction_of_tumors_for_comparison*num_tumors
dist <- 1 - as.matrix(cor_weighted)

colnames(dist) <- colnames(composite_mat)
rownames(dist) <- colnames(composite_mat)
dist_tumors_only <- dist[setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA),setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)]

# Calculate the standard deviations for categorization 
median_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
mad_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
mean_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
sd_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
for(i in 1:num_tumors){
  median_dist_tumor_to_k_nearest_tumors[i] <- median(sort(dist_tumors_only[i,-i])[1:k])
  mad_dist_tumor_to_k_nearest_tumors[i] <- mad(sort(dist_tumors_only[i,-i])[1:k])
  mean_dist_tumor_to_k_nearest_tumors[i] <- mean(sort(dist_tumors_only[i,-i])[1:k])
  sd_dist_tumor_to_k_nearest_tumors[i] <- sd(sort(dist_tumors_only[i,-i])[1:k])
}

dist_cell_line_to_nearest_tumor <- rep(NA, num_cell_lines)
median_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
mad_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
mean_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
sd_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
dist_cell_line_to_nearest_tumor <- rep(NA, num_cell_lines)
for(i in 1:num_cell_lines) {
  mad_dist_cell_line_to_k_nearest_tumors[i] <- mad(sort(dist[cell_lines_with_both_MUT_and_CNA[i],setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)])[1:k])
  median_dist_cell_line_to_k_nearest_tumors[i] <- median(sort(dist[cell_lines_with_both_MUT_and_CNA[i],setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)])[1:k])
  sd_dist_cell_line_to_k_nearest_tumors[i] <- sd(sort(dist[cell_lines_with_both_MUT_and_CNA[i],setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)])[1:k])
  mean_dist_cell_line_to_k_nearest_tumors[i] <- mean(sort(dist[cell_lines_with_both_MUT_and_CNA[i],setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)])[1:k])
  dist_cell_line_to_nearest_tumor[i] <- min(dist[cell_lines_with_both_MUT_and_CNA[i],setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)])
}
names(median_dist_cell_line_to_k_nearest_tumors) <- cell_lines_with_both_MUT_and_CNA
names(mean_dist_cell_line_to_k_nearest_tumors) <- cell_lines_with_both_MUT_and_CNA
names(mad_dist_cell_line_to_k_nearest_tumors) <- cell_lines_with_both_MUT_and_CNA
names(sd_dist_cell_line_to_k_nearest_tumors) <- cell_lines_with_both_MUT_and_CNA
names(dist_cell_line_to_nearest_tumor) <- cell_lines_with_both_MUT_and_CNA

cutoff_high <- mean(mean_dist_tumor_to_k_nearest_tumors) + 3*sd(mean_dist_tumor_to_k_nearest_tumors)
cutoff_low <- mean(mean_dist_tumor_to_k_nearest_tumors) + sd(mean_dist_tumor_to_k_nearest_tumors)
median_similarity_cell_line_to_k_nearest_tumors <- 1 - median_dist_cell_line_to_k_nearest_tumors
mean_similarity_cell_line_to_k_nearest_tumors <- 1 - mean_dist_cell_line_to_k_nearest_tumors
median_similarity_tumor_to_k_nearest_tumors <- 1 - median_dist_tumor_to_k_nearest_tumors
mean_similarity_tumor_to_k_nearest_tumors <- 1 - mean_dist_tumor_to_k_nearest_tumors
cutoff_high_similarity <-  mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors)
cutoff_low_similarity <- mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)
great_matches <- cell_line_ids[which(mean_similarity_cell_line_to_k_nearest_tumors >= mean(mean_similarity_tumor_to_k_nearest_tumors))]
good_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < mean(mean_similarity_tumor_to_k_nearest_tumors))) ]
moderately_good_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= (mean(mean_similarity_tumor_to_k_nearest_tumors) - 2*sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors))))) ]
poor_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= (mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - 2*sd(mean_similarity_tumor_to_k_nearest_tumors))))) ]
outliers <- cell_line_ids[which(mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)))]

# MDS PLOT ----
num_cell_lines <- length(intersect(colnames(composite_mat),cell_lines_with_both_MUT_and_CNA))
num_tumors <- length(intersect(colnames(composite_mat),tumors_with_both_MUT_and_CNA))
cell_lines_and_tumors.col <- c(rep("orange",num_cell_lines),rep("blue",num_tumors))
cell_lines_and_tumors.pch <- c(rep(17,num_cell_lines),rep(20,num_tumors))
freq_alt_samplewise_mut <- apply(composite_MUT,2,compute_freq_alt)
freq_alt_samplewise_cna <- apply(composite_CNA,2,compute_freq_alt)
freq_alt_samplewise_cna_high_level_only <- apply(composite_CNA_high_level_only,2,compute_freq_alt)
plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlim=c(0,1),ylim=c(0,1),xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",main="Using all CNAs")
text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_line_ids,cex=0.6)
legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))
plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_line_ids,cex=0.6)
legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))

isomdsfit <- isoMDS(1-cor_weighted,k=2)
plot(isomdsfit$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, including low-level CNAs")
text(x=isomdsfit$points[1:num_cell_lines,1], y=isomdsfit$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)

isomdsfit_dfrm <- as.data.frame(isomdsfit$points)
scatter <- ggplot(isomdsfit_dfrm, aes(V1,V2))
scatter + geom_point(colour=cell_lines_and_tumors.col, shape=cell_lines_and_tumors.pch, size=3) + labs(x = "Coordinate 1", y = "Coordinate 2") + geom_text(aes(label=c(cell_line_ids,rep("",num_tumors)),hjust=0.5, vjust=-1))

isomdsfit2 <- isoMDS(1- cor_weighted_high_level_only,k=2)
plot(isomdsfit2$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, excluding low-level CNAs")
text(x=isomdsfit2$points[1:num_cell_lines,1], y=isomdsfit2$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)
legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))


