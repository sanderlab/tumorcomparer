#' Computes frequency of alteration (fraction of non-zeroes)
#' 
#' TODO MUT
#'
#' @param alt_mat a matrix or vector of discrete values - 0s represent unaltered, non-zeroes altered
#' 
#' @return the fraction of non-zeroes, (i.e. altered features/genes)
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
compute_freq_alt <- function(alt_mat) {
  return(length(which(alt_mat != 0))/length(alt_mat)) # fraction of non-zeroes, ie altered features/genes
}