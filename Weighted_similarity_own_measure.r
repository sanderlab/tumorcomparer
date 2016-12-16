library(cluster)
library(MASS)
tumor_MUT <- read.table("tumor_MUT.txt",sep="\t",header=T,row.names=1)
tumor_CNA <- read.table("tumor_CNA.txt",sep="\t",header=T,row.names=1)
cell_line_MUT <- read.table("cell_line_MUT.txt",sep="\t",header=T,row.names=1)
cell_line_CNA <- read.table("cell_line_CNA.txt",sep="\t",header=T,row.names=1)

CNA_default_weight <- 0.01
MUT_default_weight <- 0.01

tumors_with_both_MUT_and_CNA <- intersect(colnames(tumor_MUT),colnames(tumor_CNA))
cell_lines_with_both_MUT_and_CNA <- intersect(colnames(cell_line_MUT),colnames(cell_line_CNA))
genes_with_MUT_in_both <- intersect(rownames(tumor_MUT),rownames(cell_line_MUT))
genes_with_CNA_in_both <- intersect(rownames(tumor_CNA),rownames(cell_line_CNA))
genes_in_all_4_files <- intersect(genes_with_MUT_in_both,genes_with_CNA_in_both) # Need not do this unless really want MUT and CNA data for the same gene sets

keep_only_high_level_cnas <- function(cna_mat)
{
  cna_mat_high_only <- cna_mat
  cna_mat_high_only[which(abs(cna_mat) == 1)] <- 0
  return(cna_mat_high_only)
}

cell_line_CNA_high_level_only <- apply(cell_line_CNA,2,keep_only_high_level_cnas)
tumor_CNA_high_level_only <- apply(tumor_CNA,2,keep_only_high_level_cnas)

#composite_MUT <- cbind(cell_line_MUT[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_MUT[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
composite_CNA <- cbind(cell_line_CNA[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
composite_CNA_high_level_only <- cbind(cell_line_CNA_high_level_only[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA_high_level_only[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
composite_MUT <- cbind(cell_line_MUT[genes_with_MUT_in_both,cell_lines_with_both_MUT_and_CNA],tumor_MUT[genes_with_MUT_in_both,tumors_with_both_MUT_and_CNA])
#composite_CNA <- cbind(cell_line_CNA[genes_with_CNA_in_both,cell_lines_with_both_MUT_and_CNA],tumor_CNA[genes_with_CNA_in_both,tumors_with_both_MUT_and_CNA])
#composite_CNA_high_level_only <- cbind(cell_line_CNA_high_level_only[genes_with_CNA_in_both,cell_lines_with_both_MUT_and_CNA],tumor_CNA_high_level_only[genes_with_CNA_in_both,tumors_with_both_MUT_and_CNA])

rownames(composite_MUT) <- paste(rownames(composite_MUT),"MUT",sep="_")
rownames(composite_CNA) <- paste(rownames(composite_CNA),"CNA",sep="_")
rownames(composite_CNA_high_level_only) <- paste(rownames(composite_CNA_high_level_only),"CNA",sep="_")

composite_mat <- rbind(composite_MUT,composite_CNA)
alt_mat <- rbind(composite_MUT,composite_CNA)
write.table(composite_mat,file="composite_alteration_matrix.txt",sep="\t",quote=F)
composite_mat_high_level_only <- rbind(composite_MUT,composite_CNA_high_level_only)

#alt_mat <- read.table("alteration_matrix.txt",sep="\t",header=T,row.names=1)
# Assign frequency weights as (freq. of alteration of gene)/(mean freq. of alteration across all genes) - "rewarding recurrent changes"
overall_alt_freq <- length(which((alt_mat[]) != 0))/ ( length(which((alt_mat[])==0)) + length(which((alt_mat[])!=0)))

compute_freq_alt <- function(alt_mat)
{
 return(length(which(alt_mat!=0))/length(alt_mat))
}

freq_alt <- rep(0,nrow(alt_mat))
freq_alt <- apply(alt_mat,1,compute_freq_alt)

freq_alt_samplewise <- apply(alt_mat,2,compute_freq_alt)
freq_alt_samplewise_CNA_high_level_only <- apply(composite_mat_high_level_only,2,compute_freq_alt)
#alt_mat <- alt_mat[,which(freq_alt_samplewise > 0)]
#composite_mat <- composite_mat[,which(freq_alt_samplewise > 0)]
#composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise > 0)]

#alt_mat <- alt_mat[which(freq_alt > 0),] # remove genes/events in which no sample shows mutation/alteration
#freq_alt <- freq_alt[which(freq_alt > 0)]

names(freq_alt) <- rownames(alt_mat)
freq_weights <- rep(1,nrow(alt_mat))
freq_weights <- freq_alt/overall_alt_freq
names(freq_weights) <- rownames(alt_mat)

returnfirstpart <- function(x)
{
  return(strsplit(x,"_")[[1]][1])
}

cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, returnfirstpart)



# Read in user-provided weights
known_cancer_genes_and_weights_all <- read.table("Default_weights_for_known_cancer_genes.txt",sep="\t",header=T,row.names=1) 
known_cancer_genes_and_weights <- as.matrix(known_cancer_genes_and_weights_all[intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
genes_and_weights_all <- read.table("Genes_and_weights.txt",sep="\t",header=T,row.names=1) 
genes_and_weights <- as.matrix(genes_and_weights_all[intersect(rownames(genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
rownames(genes_and_weights) <- intersect(rownames(genes_and_weights_all),rownames(alt_mat))
rownames(known_cancer_genes_and_weights) <- intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat))
#annotation_weights <- rep(1,nrow(alt_mat))
#annotation_weights <- rep(0,nrow(alt_mat)) #Exclude all but specified genes
annotation_weights <- rep(0.1,nrow(alt_mat)) 
#annotation_weights <- freq_weights # initialize with relative freq. of alterations as weight
#annotation_weights[grep("_CNA",rownames(alt_mat))] <- freq_weights[grep("_CNA",rownames(alt_mat))]
annotation_weights[grep("_CNA",rownames(alt_mat))] <- CNA_default_weight
annotation_weights[grep("_MUT",rownames(alt_mat))] <- MUT_default_weight
#min_wt <- 1/nrow(alt_mat)
#annotation_weights <- rep(min_wt,nrow(alt_mat)) 
names(annotation_weights) <- rownames(alt_mat) 
for(i in 1:nrow(known_cancer_genes_and_weights)) # Overwrite default weights for known cancer genes
    annotation_weights[rownames(known_cancer_genes_and_weights)[i]] = known_cancer_genes_and_weights[i,]
for(i in 1:nrow(genes_and_weights)) # Overwrite default weights for input provided weights
    annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i,]


gene_weights <- rep(1,nrow(alt_mat))
names(gene_weights) <- rownames(alt_mat)
#gene_weights <- freq_weights * annotation_weights # combine freq. weights and user-provided weights
gene_weights <- annotation_weights # if using user-provided weights only

#alt_mat <- alt_mat[which(freq_alt > 0),] # remove genes/events in which no sample shows mutation/alteration
#freq_alt <- freq_alt[which(freq_alt > 0)]
#gene_weights <- gene_weights[which(freq_alt > 0)]

#gene_weights <- gene_weights/sum(gene_weights) # normalize weights
gene_weights <- gene_weights/max(gene_weights) # map to 0-1

pseudoweight <- 0.000001
#standardize_min_max <- function(x){return( (x-min(x))/(max(x)-min(x))  )} # standardize to [0,1] by subtracting minimum and dividing by range
#standardize_min_max <- function(x) # standardize to [0,1] by subtracting minimum and dividing by range
#{
#  tmp <-  (x-min(x)+pseudoweight)/(max(x)-min(x)+pseudoweight)  
#  return(tmp-min(tmp))
#}
#alt_mat_std <- t(apply((as.matrix(alt_mat)),1,standardize_min_max))
#alt_mat <- alt_mat_std


weighted.corr <- function (a, b, w = rep(1, nrow(a))/nrow(a))
{
    # normalize weights
    w <- w / sum(w)

    # center matrices
    a <- sweep(a, 2, colSums(a * w))
    b <- sweep(b, 2, colSums(b * w))

    # compute weighted correlation
    t(w*a) %*% b / sqrt( colSums(w * a**2) %*% t(colSums(w * b**2)) )
}

cor_weighted <- weighted.corr(as.matrix(alt_mat),as.matrix(alt_mat),gene_weights)
cor_unweighted <- cor(alt_mat)

num_cell_lines <- length(cell_lines_with_both_MUT_and_CNA)
num_tumors <- length(tumors_with_both_MUT_and_CNA)
cell_lines_and_tumors.col <- c(rep("orange",num_cell_lines),rep("blue",num_tumors))
cell_lines_and_tumors.pch <- c(rep(17,num_cell_lines),rep(20,num_tumors))
freq_alt_samplewise_mut <- apply(composite_MUT,2,compute_freq_alt)
freq_alt_samplewise_cna <- apply(composite_CNA,2,compute_freq_alt)
freq_alt_samplewise_cna_high_level_only <- apply(composite_CNA_high_level_only,2,compute_freq_alt)
plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlim=c(0,1),ylim=c(0,1),xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",main="Using all CNAs")
text(freq_alt_samplewise_mut[1:length(cell_line_ids)],freq_alt_samplewise_cna[1:length(cell_line_ids)],labels=cell_line_ids,cex=0.6)
plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_lines_with_both_MUT_and_CNA,cex=0.6)
#text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_line_ids,cex=0.6)

compute_weighted_distance_excluding_zero_zero_matches <- function(alt_mat, gene_weights)
{
  for(i in 1:ncol(alt_mat))
  for(j in 1:ncol(alt_mat))
  {  
    if(sum(gene_weights[which(((alt_mat[,i]!=0) | (alt_mat[,j]!=0)))]) == 0) 
    {
      weighted_distance_excluding_zero_zero_matches[i,j] <- 1e-6 # or some default
    }
    else
    {
      weighted_distance_excluding_zero_zero_matches[i,j] <- 1e-6 + sum(gene_weights[which(alt_mat[,i] != alt_mat[,j])])/sum(gene_weights[which(((alt_mat[,i]!=0) | (alt_mat[,j]!=0)))])
    }
  } 
    return(weighted_distance_excluding_zero_zero_matches)
}

compute_weighted_distance_including_zero_zero_matches <- function(alt_mat, gene_weights)
{
  for(i in 1:ncol(alt_mat))
  for(j in 1:ncol(alt_mat))
  {
      weighted_distance_including_zero_zero_matches[i,j] <- sum(gene_weights[which(alt_mat[,i] != alt_mat[,j])])/sum(gene_weights)
      #weighted_distance_including_zero_zero_matches[i,j] <- length(which(alt_mat[,i] != alt_mat[,j]))/nrow(alt_mat)
  }
  return(weighted_distance_including_zero_zero_matches)
}



weighted_distance_excluding_zero_zero_matches <- matrix(NA,nrow=ncol(alt_mat),ncol=ncol(alt_mat))
weighted_distance_excluding_zero_zero_matches_MUT_only <- matrix(NA,nrow=ncol(alt_mat),ncol=ncol(alt_mat))
weighted_distance_excluding_zero_zero_matches_high_level_only <- matrix(NA,nrow=ncol(alt_mat),ncol=ncol(alt_mat))
weighted_distance_excluding_zero_zero_matches_CNA_only_high_level_only <- matrix(NA,nrow=ncol(alt_mat),ncol=ncol(alt_mat))
weighted_distance_including_zero_zero_matches <- matrix(NA,nrow=ncol(alt_mat),ncol=ncol(alt_mat))

weighted_distance_excluding_zero_zero_matches_high_level_only <- compute_weighted_distance_excluding_zero_zero_matches(composite_mat_high_level_only,gene_weights)
weighted_distance_excluding_zero_zero_matches <- compute_weighted_distance_excluding_zero_zero_matches(composite_mat,gene_weights)
#weighted_distance_excluding_zero_zero_matches_MUT_only <- compute_weighted_distance_excluding_zero_zero_matches(composite_MUT,gene_weights[grep("_MUT",rownames(composite_mat))])
#weighted_distance_excluding_zero_zero_matches_CNA_only <- compute_weighted_distance_excluding_zero_zero_matches(composite_CNA,gene_weights[grep("_CNA",rownames(composite_mat))])
#weighted_distance_excluding_zero_zero_matches_CNA_only_high_level_only <- compute_weighted_distance_excluding_zero_zero_matches(composite_CNA_high_level_only,gene_weights[grep("_CNA",rownames(composite_mat))])

#weighted_distance_excluding_zero_zero_matches_MUT_only[which(is.na(weighted_distance_excluding_zero_zero_matches_MUT_only))] <- 0
#weighted_distance_excluding_zero_zero_matches_MUT_only <- weighted_distance_excluding_zero_zero_matches_MUT_only + 1e-6
#isomdsfit <- isoMDS(weighted_distance_excluding_zero_zero_matches_MUT_only,k=2)
#plot(isomdsfit$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted distance using mutations only, excluding zero-zero matches")
#text(isomdsfit$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)

#weighted_distance_excluding_zero_zero_matches_CNA_only[which(is.na(weighted_distance_excluding_zero_zero_matches_CNA_only))] <- 0
#weighted_distance_excluding_zero_zero_matches_CNA_only <- weighted_distance_excluding_zero_zero_matches_CNA_only + 1e-6
#isomdsfit <- isoMDS(weighted_distance_excluding_zero_zero_matches_CNA_only,k=2)
#plot(isomdsfit$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted distance using CNAs only, excluding zero-zero matches")
text(isomdsfit$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)


weighted_distance_excluding_zero_zero_matches[which(is.na(weighted_distance_excluding_zero_zero_matches))] <- 0
weighted_distance_excluding_zero_zero_matches <- weighted_distance_excluding_zero_zero_matches + 1e-6
isomdsfit <- isoMDS(weighted_distance_excluding_zero_zero_matches,k=2)
plot(isomdsfit$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted distance, excluding zero-zero matches")
text(isomdsfit$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)

weighted_distance_excluding_zero_zero_matches_high_level_only[which(is.na(weighted_distance_excluding_zero_zero_matches_high_level_only))] <- 0
weighted_distance_excluding_zero_zero_matches_high_level_only <- weighted_distance_excluding_zero_zero_matches_high_level_only + 1e-6
isomdsfit_high_level_only <- isoMDS(weighted_distance_excluding_zero_zero_matches_high_level_only,k=2)
plot(isomdsfit_high_level_only$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted distance, excluding zero-zero matches, excluding low-level CNAs")
text(isomdsfit_high_level_only$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)

return_top_n_shared_features <- function(alt_mat,i,j,n,gene_weights)
	return(rev(sort(gene_weights[(which((alt_mat[,i] == alt_mat[,j]) &  (alt_mat[,i] != 0) ))]))[1:n])

return_top_n_features_of_sample <- function(alt_mat,i,n,gene_weights)
        return(rev(sort(gene_weights[(which(alt_mat[,i] != 0))]))[1:n])
