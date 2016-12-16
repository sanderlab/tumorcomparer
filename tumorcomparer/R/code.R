#' Run a comparison between 
get_tumor_comparison <- function(x) {
  
}

#' 
return_first_part <- function(x) {
  return(strsplit(x,"_")[[1]][1])
}

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

# MAPPING TO COLORS (ON ORANGE-BLUE SCALE)
map_msk_to_colors <- function(mean_similarity_cell_line_to_k_nearest_tumors,mean_similarity_tumor_to_k_nearest_tumors,col1,col2,numshades) {
  maxsim <- max(mean_similarity_tumor_to_k_nearest_tumors)
  #numshades <- length(mean_similarity_tumor_to_k_nearest_tumors)
  #pal2 <- colorRampPalette(c("red","green"))
  #numshades <- 100
  #pal2 <- colorRampPalette(c("orange","blue"))
  pal2 <- colorRampPalette(c(col1,col2))
  colors2 <- pal2(numshades)
  morecolors <- rep(NA,length(mean_similarity_cell_line_to_k_nearest_tumors))
  for(i in 1:length(morecolors))
  {
    colindex <- round((mean_similarity_cell_line_to_k_nearest_tumors[i]/maxsim)*numshades)
    if(colindex==0){colindex <- 1}
    if(colindex > numshades){colindex <- numshades}
    morecolors[i] <- colors2[colindex]
  }  
  return(morecolors)
}
