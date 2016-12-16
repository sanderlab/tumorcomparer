#' X 
keep_only_high_level_cnas <- function(cna_mat) {
  cna_mat_high_only <- cna_mat
  cna_mat_high_only[which(abs(cna_mat) == 1)] <- 0
  return(cna_mat_high_only)
}