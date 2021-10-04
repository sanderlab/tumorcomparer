#' Calculate a weighted distance (between two vectors) based on Jaccard's coefficient
#' 
#' @param x a vector of discrete values
#' @param y a vector of discrete values
#' @param weights a vector of non-negative alteration values
#' 
#' @return a weighted distance based on Jaccard's coefficient. Ratio of the size 
#'   of intersection to the size of the union of two sets after discarding 0-0 
#'   matches from the two vectors
#'   
#' @examples 
#' n <- 100
#' x <- sample(c(0, 1), n, replace=TRUE)
#' y <- sample(c(0, 1), n, replace=TRUE)
#' weights <- rnorm(n)
#' pair_dist(x, y, weights)
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
pair_dist <- function(x, y, weights) {
  # sum of weights of alterations shared by both samples, divided by sum of
  # weights of alterations present in at least one of the samples
  (1e-6 + sum(weights[which(x != y)]))/(1e-6 + sum(weights[which(((x!=0) | (y!=0)))]))
}
