#' Map mean similarity of cell lines to k nearest tumors to create a gradient from color 1 to color 2
#'
#' @param mean_similarity_cell_line_to_k_nearest_tumors a vector of mean similarities for each cell line to k nearest tumors
#' @param mean_similarity_tumor_to_k_nearest_tumors a vector of mean similarities for each tumor to k nearest tumors
#' @param col1 a color (DEFAULT: orange)
#' @param col2 a color (DEFAULT: blue)
#' @param numshades the number of shades on the color scale (DEFAULT: 100)
#' 
#' @return a vector of colors of length equal to mean_similarity_cell_line_to_k_nearest_tumors
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
map_mean_similarity_to_gradient <- function(mean_similarity_cell_line_to_k_nearest_tumors,
           mean_similarity_tumor_to_k_nearest_tumors,
           col1 = "orange",
           col2 = "blue",
           numshades = 100) {
    
    maxsim <- max(mean_similarity_tumor_to_k_nearest_tumors)
    
    #numshades <- length(mean_similarity_tumor_to_k_nearest_tumors)
    #pal2 <- colorRampPalette(c("red","green"))
    #numshades <- 100
    #pal2 <- colorRampPalette(c("orange","blue"))
    
    pal2 <- colorRampPalette(c(col1, col2))
    colors2 <- pal2(numshades)
    morecolors <- rep(NA, length(mean_similarity_cell_line_to_k_nearest_tumors))
    
    # Get the right number of colors 
    for (i in 1:length(morecolors)) {
      colindex <- round((mean_similarity_cell_line_to_k_nearest_tumors[i] / maxsim) * numshades)
      if (colindex == 0) {
        colindex <- 1
      }
      if (colindex > numshades) {
        colindex <- numshades
      }
      morecolors[i] <- colors2[colindex]
    }
    
    return(morecolors)
  }
