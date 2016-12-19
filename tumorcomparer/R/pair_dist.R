#' Calculate pairwise distance [FIX]
#' 
#' @param x a matrix [FIX]
#' @param y a matrix [FIX]
#' 
#' @return weights [FIX]
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
pair_dist <- function(x, y, weights) {
  (1e-6 + sum(weights[which(x != y)]))/(1e-6 + sum(weights[which(((x!=0) | (y!=0)))]))
}