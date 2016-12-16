#' 
calc_weighted_corr <- function (a, b, w = rep(1, nrow(a))/nrow(a)) {
  # normalize weights
  w <- w / sum(w)
  
  # center matrices
  a <- sweep(a, 2, colSums(a * w))
  b <- sweep(b, 2, colSums(b * w))
  
  # compute weighted correlation
  t(w*a) %*% b / sqrt( colSums(w * a**2) %*% t(colSums(w * b**2)) )
}