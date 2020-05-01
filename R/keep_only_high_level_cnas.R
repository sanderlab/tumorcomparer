#' Keep only high-level CNAs (set low-level CNAs to 0)
#' 
#' @param cna_mat a matrix or vector of 5-valued CNA (copy number alteration) data; 
#'   possible values: -2: Deep deletion, -1: Shallow deletion, 0: Diploid (or default copy number), 
#'   1: Low-level gain, 2: High-level amplification
#'   
#' @return a matrix where low-level gains have been set to 0
#' 
#' @examples 
#' set.seed(1)
#' tmp <- sample(c(-2, -1, 0, 1, 2), 100, replace=TRUE)
#' mat <- matrix(tmp, 10, 10)
#' keep_only_high_level_cnas(mat)
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
keep_only_high_level_cnas <- function(cna_mat) {
  cna_mat_high_only <- cna_mat # copy input matrix
  cna_mat_high_only[which(abs(cna_mat) == 1)] <- 0 # set low-level CNAs in copy to 0
  
  return(cna_mat_high_only) # return copy
}