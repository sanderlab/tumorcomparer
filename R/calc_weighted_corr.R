#' Calculate weighted correlation 
#' 
#' @param a numeric matrix
#' @param b numeric matrix 
#' @param w alteration weights
#' 
#' @return a weighted correlation (similarity) matrix 
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @examples 
#' t1 <- c(1,1,1,1,1,1,2,2,2,3,3,3,4,4)
#' t2 <- rev(t1)
#' weight <- c(.5,.5,.5,.5,.5,1,1,1,1,2,2,2,2,2)
#' calc_weighted_corr(t1, t2, weight)
#'
#' @concept tumorcomparer
#' @export
#'
#' @importFrom weights wtd.cors 
calc_weighted_corr <- function (a, b, w = rep(1, nrow(a))/nrow(a)) {
  # normalize weights
  #w <- w / sum(w)
  
  # center matrices
  #a <- sweep(a, 2, colSums(a * w))
  #b <- sweep(b, 2, colSums(b * w))
  
  # compute weighted correlation
  #t(w*a) %*% b / sqrt( colSums(w * a**2) %*% t(colSums(w * b**2)))
  wtd.cors(x=a, y=b,  weight=w)
}
