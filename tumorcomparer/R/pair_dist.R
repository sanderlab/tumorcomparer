#' Calculate a weighted distance (between two vectors) based on Jaccard's coefficient
#' 
#' @param x a vector
#' @param y a vector
#' @param weights alteration weights
#' 
#' @return a weighted distance based on Jaccard's coefficient (ratio of the size of intersection to the size of union of two sets)
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
pair_dist <- function(x, y, weights) {
  (1e-6 + sum(weights[which(x != y)]))/(1e-6 + sum(weights[which(((x!=0) | (y!=0)))])) # sum of weights of alterations shared by both samples, divided by sum of weights of alterations present in at least one of the samples
}