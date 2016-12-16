compute_freq_alt <- function(alt_mat) {
  return(length(which(alt_mat!=0))/length(alt_mat))
}