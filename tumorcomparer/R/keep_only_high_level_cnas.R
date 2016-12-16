#' Keep only high-level CNAs (set low-level CNAs to 0)
#' @cna_mat a vector of 5-valued CNA (copy number alteration) data - 
#' -2: Deep deletion
#' -1: Shallow deletion
#' 0: Diploid (or default copy number)
#' 1: Low-level gain
#' 2: High-level amplification
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (augustin@mail.nih.gov)
#'
#' 
#'
#' @concept tumorcomparer
#' @export
keep_only_high_level_cnas <- function(cna_mat) {
  cna_mat_high_only <- cna_mat # copy input matrix
  cna_mat_high_only[which(abs(cna_mat) == 1)] <- 0 # set low-level CNAs in copy to 0
  return(cna_mat_high_only) # return copy
}