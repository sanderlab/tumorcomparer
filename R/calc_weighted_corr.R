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
#' @concept tumorcomparer
#' @export
#' @importFrom expss w_cor w_spearman w_pearson


calc_weighted_corr <- function (a, b, w = rep(1, nrow(a))/nrow(a)) {
  # normalize weights
  w <- w / sum(w)
  
  # center matrices
  a <- sweep(a, 2, colSums(a * w))
  b <- sweep(b, 2, colSums(b * w))
  
  # compute weighted correlation
  t(w*a) %*% b / sqrt( colSums(w * a**2) %*% t(colSums(w * b**2)))
  #cov.wt(cbind(a,b), wt = w, cor = TRUE)
  #w_spearman(cbind(a, b), weight = w, use = "pairwise.complete.obs")
}
