# PARAM
fraction_of_tumors_for_comparison <- 0.1 #X

cor_weighted # MAT
composite_mat # MAT
cell_lines_with_both_MUT_and_CNA # VECTOR
tumors_with_both_MUT_and_CNA #VECT

categorize_cell_lines <- function(x) {
  # CATEGORIZE ----
  num_cell_lines <- length(intersect(colnames(composite_mat),cell_lines_with_both_MUT_and_CNA))
  num_tumors <- length(intersect(colnames(composite_mat),tumors_with_both_MUT_and_CNA))
  
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
  
  results <- list(great_matches=great_matches, good_matches=good_matches, moderately_good_matches=moderately_good_matches, poor_matches=poor_matches, outliers=outliers)
  
}

