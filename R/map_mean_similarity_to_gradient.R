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
#' @examples 
#' # Generated using: tumorcomparer::run_comparison() 
#' comparison_result <- readRDS(system.file("test_output", "ov_comparison_result.rds", 
#'   package="tumorcomparer"))
#' 
#' categorization_list <- categorize_cell_lines(
#'   num_tumors_for_comparison=length(comparison_result$tumor_ids)-1, 
#'   dist_mat=comparison_result$dist_mat,
#'   cell_line_ids=comparison_result$cell_line_ids,
#'   tumor_ids=comparison_result$tumor_ids,
#'   trim_cell_line_names=FALSE) 
#'   
#' result <- map_mean_similarity_to_gradient(
#'   mean_similarity_cell_line_to_k_nearest_tumors=
#'     categorization_list$mean_similarity_cell_line_to_k_nearest_tumors,
#'   mean_similarity_tumor_to_k_nearest_tumors=
#'     categorization_list$mean_similarity_tumor_to_k_nearest_tumors,
#'   col1="orange",
#'   col2="blue", 
#'   numshades=100)
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom grDevices colorRampPalette
map_mean_similarity_to_gradient <- function(
           mean_similarity_cell_line_to_k_nearest_tumors,
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
