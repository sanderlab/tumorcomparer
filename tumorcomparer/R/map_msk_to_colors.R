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
