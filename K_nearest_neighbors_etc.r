#k <- num_tumors-1 # Number of nearest neighbors to consider
#k <- min(10, 0.5*num_tumors) # Number of nearest neighbors to consider
#k <- max(10, 0.5*num_tumors) # Number of nearest neighbors to consider
#k <- 30
k <- 0.1*num_tumors
#k <- 0.3*num_tumors
#k <- 0.9*num_tumors
#dist <- as.matrix(gowdist)
#dist <- 1 - as.matrix(cor_unweighted)
#dist <- 1 - as.matrix(cor_weighted)
#dist <- 1 - as.matrix(cosine_weighted)
#dist <- 1 - as.matrix(cor_weighted_high_level_only)
#dist <- as.matrix(unweighted_distance_excluding_zero_zero_matches)
dist <- as.matrix(weighted_distance_excluding_zero_zero_matches)
#dist <- as.matrix(weighted_distance_excluding_zero_zero_matches_high_level_only)
colnames(dist) <- colnames(composite_mat)
rownames(dist) <- colnames(composite_mat)
dist_tumors_only <- dist[setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA),setdiff(colnames(dist),cell_lines_with_both_MUT_and_CNA)]

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
