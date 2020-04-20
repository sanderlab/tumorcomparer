#' Categorize cell lines by the level of similarity to k-nearest tumors (DEFAULT: 10 percent of tumors)
#' 
#' @param num_tumors_for_comparison number of tumors used in a k-nearest
#'  neighbor comparison (DEFAULT: 10) 
#' @param dist_mat a matrix of pairwise weighted distances between all cell lines and tumors
#' @param cell_line_ids IDs/names of cell lines 
#' @param tumor_ids IDs of tumors 
#' @param trim_cell_line_names a boolean whether to trim the the cell lines; 
#'   this is optional and used for long cell line names in CCLE format (i.e. CELLLINE_TISSUE); (DEFAULT: FALSE)
#' 
#' @return a list with the following items: 
#' \itemize{
#'   \item{"mean_similarity_cell_line_to_k_nearest_tumors"}{the mean similarity of each cell line to the k-nearest tumors calculated using k=fraction_of_tumors_for_comparison*number of tumors}
#'   \item{"mean_similarity_tumor_to_k_nearest_tumors"}{the mean similarity of each tumor sample to the k-nearest tumors calculated using k=fraction_of_tumors_for_comparison*number of tumors}
#'   \item{"categorization"}{a 2-column data.frame with the cell line categorizations: Sample_ID and Category; Category values can be: "Great", "Good", "Moderately Good", "Poor", "Outliers"}
#' }
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom stats cor dist mad median sd
categorize_cell_lines <- function(num_tumors_for_comparison= 10, 
                                  dist_mat,
                                  cell_line_ids, 
                                  tumor_ids,
                                  trim_cell_line_names=FALSE) {
  
#categorize_cell_lines <- function(fraction_of_tumors_for_comparison=0.1, 
#                                  dist_mat,
#                                  composite_mat, 
#                                  cell_line_ids, 
#                                  tumor_ids,
#                                  trim_cell_line_names=FALSE) {
  
  num_cell_lines <- length(cell_line_ids)
  num_tumors <- length(tumor_ids)
  
  ## Set K-nearest neighbors 
  k <- num_tumors_for_comparison
  
  #colnames(dist_mat) <- colnames(composite_mat)
  #rownames(dist_mat) <- colnames(composite_mat)
  dist_tumors_only <- dist_mat[setdiff(colnames(dist_mat),cell_line_ids),setdiff(colnames(dist_mat),cell_line_ids)]
  
  # Calculate the standard deviations for categorization 
  median_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
  mad_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
  mean_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
  sd_dist_tumor_to_k_nearest_tumors <- rep(NA, num_tumors)
  for(i in 1:num_tumors){
    # Sort the distances of this tumor to other tumors, and save distances to the k-nearest tumors (in variable tmp)
    tmp <- sort(dist_tumors_only[i,-i])[1:k]
    
    median_dist_tumor_to_k_nearest_tumors[i] <- median(tmp)
    mad_dist_tumor_to_k_nearest_tumors[i] <- mad(tmp)
    mean_dist_tumor_to_k_nearest_tumors[i] <- mean(tmp)
    sd_dist_tumor_to_k_nearest_tumors[i] <- sd(tmp)
  }
  
  dist_cell_line_to_nearest_tumor <- rep(NA, num_cell_lines)
  median_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
  mad_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
  mean_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
  sd_dist_cell_line_to_k_nearest_tumors <- rep(NA, num_cell_lines)
  dist_cell_line_to_nearest_tumor <- rep(NA, num_cell_lines)
  for(i in 1:num_cell_lines) {
    # Sort the distances of this cell line to tumors, and save distances to the k-nearest tumors (in variable tmp)
    tmp <- sort(dist_mat[cell_line_ids[i], setdiff(colnames(dist_mat),cell_line_ids)])[1:k]

    mad_dist_cell_line_to_k_nearest_tumors[i] <- mad(tmp)
    median_dist_cell_line_to_k_nearest_tumors[i] <- median(tmp)
    sd_dist_cell_line_to_k_nearest_tumors[i] <- sd(tmp)
    mean_dist_cell_line_to_k_nearest_tumors[i] <- mean(tmp)
    
    # Compute distance to nearest tumor
    #dist_cell_line_to_nearest_tumor[i] <- min(dist_mat[cell_line_ids[i],setdiff(colnames(dist_mat),cell_line_ids)])
  }
  
  names(median_dist_cell_line_to_k_nearest_tumors) <- cell_line_ids
  names(mean_dist_cell_line_to_k_nearest_tumors) <- cell_line_ids
  names(mad_dist_cell_line_to_k_nearest_tumors) <- cell_line_ids
  names(sd_dist_cell_line_to_k_nearest_tumors) <- cell_line_ids
  names(dist_cell_line_to_nearest_tumor) <- cell_line_ids
  
  #cutoff_high <- mean(mean_dist_tumor_to_k_nearest_tumors) + 3*sd(mean_dist_tumor_to_k_nearest_tumors)
  #cutoff_low <- mean(mean_dist_tumor_to_k_nearest_tumors) + sd(mean_dist_tumor_to_k_nearest_tumors)
  
  median_similarity_cell_line_to_k_nearest_tumors <- 1 - median_dist_cell_line_to_k_nearest_tumors
  mean_similarity_cell_line_to_k_nearest_tumors <- 1 - mean_dist_cell_line_to_k_nearest_tumors
  median_similarity_tumor_to_k_nearest_tumors <- 1 - median_dist_tumor_to_k_nearest_tumors
  mean_similarity_tumor_to_k_nearest_tumors <- 1 - mean_dist_tumor_to_k_nearest_tumors
  #cutoff_high_similarity <-  mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors)
  #cutoff_low_similarity <- mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)
  
  if(trim_cell_line_names) {
    cell_line_ids <- sapply(cell_line_ids, return_first_part)
  } 
  
  # Categorize matches 
  great_matches <- cell_line_ids[which(mean_similarity_cell_line_to_k_nearest_tumors >= mean(mean_similarity_tumor_to_k_nearest_tumors))]
  good_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < mean(mean_similarity_tumor_to_k_nearest_tumors))) ]
  moderately_good_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= (mean(mean_similarity_tumor_to_k_nearest_tumors) - 2*sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - sd(mean_similarity_tumor_to_k_nearest_tumors))))) ]
  poor_matches <- cell_line_ids[which((mean_similarity_cell_line_to_k_nearest_tumors >= (mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)) & (mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - 2*sd(mean_similarity_tumor_to_k_nearest_tumors))))) ]
  outliers <- cell_line_ids[which(mean_similarity_cell_line_to_k_nearest_tumors < (mean(mean_similarity_tumor_to_k_nearest_tumors) - 3*sd(mean_similarity_tumor_to_k_nearest_tumors)))]

  matches_samples <- c(great_matches, good_matches, moderately_good_matches, poor_matches, outliers)
  matches_categories <- c(rep("Great", length(great_matches)), 
                      rep("Good", length(good_matches)), 
                      rep("Moderately Good", length(moderately_good_matches)), 
                      rep("Poor", length(poor_matches)),
                      rep("Outlier", length(outliers))
  )
  
  # Create single data.frame
  categorization <- data.frame(Sample_ID= matches_samples, Category=matches_categories, stringsAsFactors = FALSE)
  
  # Create results
  results <- list(
    mean_similarity_cell_line_to_k_nearest_tumors=mean_similarity_cell_line_to_k_nearest_tumors, 
    mean_similarity_tumor_to_k_nearest_tumors=mean_similarity_tumor_to_k_nearest_tumors,
    categorization=categorization
  )

  return(results)
}
