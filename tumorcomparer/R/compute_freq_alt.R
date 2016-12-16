#'Computes frequency of alteration (fraction of non-zeroes)
#' @alt_mat a vector of discrete values - 0s represent unaltered, non-zeroes altered
#' 
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (augustin@mail.nih.gov)
#'
#' 
#'
#' @concept tumorcomparer
#' @export
compute_freq_alt <- function(alt_mat) {
  return(length(which(alt_mat!=0))/length(alt_mat)) # fraction of non-zeroes, ie altered features/genes
}