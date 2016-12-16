library(cluster)
library(MASS)
library(ggplot2)

# LOAD DATA ---- 
tumor_MUT <- read.table("tumor_MUT.txt",sep="\t",header=T,row.names=1)
tumor_CNA <- read.table("tumor_CNA.txt",sep="\t",header=T,row.names=1)
cell_line_MUT <- read.table("cell_line_MUT.txt",sep="\t",header=T,row.names=1)
cell_line_CNA <- read.table("cell_line_CNA.txt",sep="\t",header=T,row.names=1)

tumors_with_both_MUT_and_CNA <- intersect(colnames(tumor_MUT),colnames(tumor_CNA))
cell_lines_with_both_MUT_and_CNA <- intersect(colnames(cell_line_MUT),colnames(cell_line_CNA))
genes_with_MUT_in_both <- intersect(rownames(tumor_MUT),rownames(cell_line_MUT))
genes_with_CNA_in_both <- intersect(rownames(tumor_CNA),rownames(cell_line_CNA))
genes_in_all_4_files <- intersect(genes_with_MUT_in_both,genes_with_CNA_in_both) # Need not do this unless really want MUT and CNA data for the same gene sets

CNA_default_weight <- 0.01
MUT_default_weight <- 0.01

CNA_known_cancer_gene_weight <- 0.1
MUT_known_cancer_gene_weight <- 0.1

keep_only_high_level_cnas <- function(cna_mat)
{
  cna_mat_high_only <- cna_mat
  cna_mat_high_only[which(abs(cna_mat) == 1)] <- 0
  return(cna_mat_high_only)
}

cell_line_CNA_high_level_only <- apply(cell_line_CNA,2,keep_only_high_level_cnas)
tumor_CNA_high_level_only <- apply(tumor_CNA,2,keep_only_high_level_cnas)

# CREATE COMPOSITE MATRIX ----
## Rows genes (MUT, CNA), columns (Tumors/cell lines)
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
freq_alt_high_level  <- apply(composite_mat_high_level_only,1,compute_freq_alt)
freq_alt_mut_tumors <- apply(composite_MUT[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)
freq_alt_cna_tumors <- apply(composite_CNA[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)

freq_alt_samplewise <- apply(alt_mat,2,compute_freq_alt)
freq_alt_samplewise_CNA_high_level_only <- apply(composite_mat_high_level_only,2,compute_freq_alt)
alt_mat <- alt_mat[,which(freq_alt_samplewise > 0)]
composite_mat <- composite_mat[,which(freq_alt_samplewise > 0)]
composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise > 0)]
#composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise_CNA_high_level_only > 0)]

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

# GET WEIGHTS 
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
#annotation_weights[grep("_CNA",rownames(alt_mat))] <- freq_alt_cna_tumors
#annotation_weights[grep("_MUT",rownames(alt_mat))] <- freq_alt_mut_tumors
#min_wt <- 1/nrow(alt_mat)
#annotation_weights <- rep(min_wt,nrow(alt_mat)) 
names(annotation_weights) <- rownames(alt_mat) 
for(i in 1:nrow(known_cancer_genes_and_weights)) # Overwrite default weights for known cancer genes
    annotation_weights[rownames(known_cancer_genes_and_weights)[i]] = known_cancer_genes_and_weights[i,]
for(i in 1:nrow(genes_and_weights)) # Overwrite default weights for input provided weights
#    annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i,] + annotation_weights[rownames(genes_and_weights)[i]]
    annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i,]


gene_weights <- rep(1,nrow(alt_mat))
names(gene_weights) <- rownames(alt_mat)
#gene_weights <- freq_weights * annotation_weights # combine freq. weights and user-provided weights
gene_weights <- annotation_weights # if using user-provided weights only

#alt_mat <- alt_mat[which(freq_alt > 0),] # remove genes/events in which no sample shows mutation/alteration
#composite_mat <- composite_mat[which(freq_alt > 0),] # remove genes/events in which no sample shows mutation/alteration
#composite_mat_high_level_only <- composite_mat_high_level_only[which(freq_alt_high_level > 0),] # remove genes/events in which no sample shows mutation/alteration
#gene_weights <- gene_weights[which(freq_alt > 0)]
#gene_weights_high_level <- gene_weights[which(freq_alt_high_level > 0)]
#freq_alt <- freq_alt[which(freq_alt > 0)]
#freq_alt_high_level <- freq_alt_high_level[which(freq_alt_high_level > 0)]

#gene_weights <- gene_weights/sum(gene_weights) # normalize weights
gene_weights <- gene_weights/max(gene_weights) # map to 0-1

# DEFAULT WEIGHT
pseudoweight <- 0.000001
#standardize_min_max <- function(x){return( (x-min(x))/(max(x)-min(x))  )} # standardize to [0,1] by subtracting minimum and dividing by range
#standardize_min_max <- function(x) # standardize to [0,1] by subtracting minimum and dividing by range
#{
#  tmp <-  (x-min(x)+pseudoweight)/(max(x)-min(x)+pseudoweight)  
#  return(tmp-min(tmp))
#}
#alt_mat_std <- t(apply((as.matrix(alt_mat)),1,standardize_min_max))
#alt_mat <- alt_mat_std

# WEIGHTED-CORRELATION FUNCTION ----
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

# CALCULATE CORRELATIONS ----
# Including low-level CNAs
cor_weighted <- weighted.corr(as.matrix(composite_mat),as.matrix(composite_mat),gene_weights)
# Excluding low-levels CNAs
cor_weighted_high_level_only <- weighted.corr(as.matrix(composite_mat_high_level_only),as.matrix(composite_mat_high_level_only),gene_weights)
cor_unweighted <- cor(alt_mat)

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
#text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_lines_with_both_MUT_and_CNA,cex=0.6)
text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_line_ids,cex=0.6)
legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))


isomdsfit <- isoMDS(1-cor_weighted,k=2)
plot(isomdsfit$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, including low-level CNAs")
#text(isomdsfit$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)
text(x=isomdsfit$points[1:num_cell_lines,1], y=isomdsfit$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)
#legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))

isomdsfit_dfrm <- as.data.frame(isomdsfit$points)
scatter <- ggplot(isomdsfit_dfrm, aes(V1,V2))
scatter + geom_point(colour=cell_lines_and_tumors.col, shape=cell_lines_and_tumors.pch, size=3) + labs(x = "Coordinate 1", y = "Coordinate 2") + geom_text(aes(label=c(cell_line_ids,rep("",num_tumors)),hjust=0.5, vjust=-1))

isomdsfit2 <- isoMDS(1- cor_weighted_high_level_only,k=2)
plot(isomdsfit2$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, excluding low-level CNAs")
#text(isomdsfit$points[1:num_cell_lines,],labels=cell_line_ids,cex=0.6)
text(x=isomdsfit2$points[1:num_cell_lines,1], y=isomdsfit2$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)
legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))

# START IGNORE ----
return_top_n_shared_features <- function(alt_mat,i,j,n,gene_weights)
	return(rev(sort(gene_weights[(which((alt_mat[,i] == alt_mat[,j]) &  (alt_mat[,i] != 0) ))]))[1:n])

return_top_n_features_of_sample <- function(alt_mat,i,n,gene_weights)
        return(rev(sort(gene_weights[(which(alt_mat[,i] != 0))]))[1:n])
# END IGNORE ----

# CATEGORIZE ----
## Set K-nearest neighbors 
#k <- num_tumors-1 # Number of nearest neighbors to consider
#k <- min(10, 0.5*num_tumors) # Number of nearest neighbors to consider
#k <- max(10, 0.5*num_tumors) # Number of nearest neighbors to consider
#k <- 30
k <- 0.1*num_tumors
#k <- 0.3*num_tumors
#k <- 0.9*num_tumors
#dist <- as.matrix(gowdist)
#dist <- 1 - as.matrix(cor_unweighted)
dist <- 1 - as.matrix(cor_weighted)
#dist <- 1 - as.matrix(cosine_weighted)
#dist <- 1 - as.matrix(cor_weighted_high_level_only)
#dist <- as.matrix(unweighted_distance_excluding_zero_zero_matches)
#dist <- as.matrix(weighted_distance_excluding_zero_zero_matches)
#dist <- as.matrix(weighted_distance_excluding_zero_zero_matches_high_level_only)
colnames(dist) <- colnames(composite_mat)
rownames(dist) <- colnames(composite_mat)
dist_tumors_only <- dist[setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA),setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)]

# Calculate the standard deviations for categorization 
median_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
mad_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
mean_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
sd_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
for(i in 1:num_tumors) 
{
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
for(i in 1:num_cell_lines)
{
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

# MAPPING TO COLORS (ON ORANGE-BLUE SCALE)
#Map_MSK_to_Colors <- function(mean_similarity_cell_line_to_k_nearest_tumors,mean_similarity_tumor_to_k_nearest_tumors)
Map_MSK_to_Colors <- function(mean_similarity_cell_line_to_k_nearest_tumors,mean_similarity_tumor_to_k_nearest_tumors,col1,col2,numshades)
{
    maxsim <- max(mean_similarity_tumor_to_k_nearest_tumors)
    #numshades <- length(mean_similarity_tumor_to_k_nearest_tumors)
    #pal2 <- colorRampPalette(c("red","green"))
    #numshades <- 100
    #pal2 <- colorRampPalette(c("orange","blue"))
    pal2 <- colorRampPalette(c(col1,col2))
    colors2 <- pal2(numshades)
    morecolors <- rep(NA,length(mean_similarity_cell_line_to_k_nearest_tumors))
    for(i in 1:length(morecolors))
    {
       colindex <- round((mean_similarity_cell_line_to_k_nearest_tumors[i]/maxsim)*numshades)
       if(colindex==0){colindex <- 1}
       if(colindex > numshades){colindex <- numshades}
       morecolors[i] <- colors2[colindex]
    }  
    return(morecolors)
}
