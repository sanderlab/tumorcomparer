#' Plot the results from run_comparison as a two dimension multi-dimensional scaling (MDS) plot
#' 
#' @param comparison_list the results of run_comparison() (See: run_comparison())
#' @param categorization_list the results of categorize_cell_lines() (See: categorize_cell_lines())
#' @param tumor_color a color for tumor points (DEFAULT: blue)
#' @param cell_line_color a color for tumor points (DEFAULT: orange)
#' @param use_gradient a boolean, cell lines will appear with a gradient of 
#'   colors with those with the cell_line_color being least similar to the tumors 
#'   and those being similar to tumors will have a color closer to the tumor_color.
#' @param tumor_shape an integer for an R plot PCH symbol (DEFAULT: 17)
#' @param cell_line_shape an integer for an R plot PCH symbol (DEFAULT: 20)
#' @param trim_cell_line_names a boolean whether to trim the the cell lines; 
#'   this is optional and used for long cell line names in CCLE format (i.e. CELLLINE_TISSUE)
#' 
#' @return Nothing is returned 
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @seealso \code{\link{run_comparison}}, \code{\link{categorize_cell_lines}} 
#'
#' @concept tumorcomparer
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text theme_minimal aes_string
plot_mds <- function(comparison_list,
                     categorization_list,
                     trim_cell_line_names=FALSE,
                     tumor_color="blue", 
                     cell_line_color="orange", 
                     use_gradient=TRUE,
                     tumor_shape=20, 
                     cell_line_shape=17) {
  # NOTE: Dist should be 
  #dist_mat <- 1-cor_weighted
  
  # From run_comparison
  isomdsfit <- comparison_list$isomdsfit
  cell_line_ids <- comparison_list$cell_line_ids
  tumor_ids <- comparison_list$tumor_ids
  
  # From categorize_cell_lines
  mean_similarity_cell_line_to_k_nearest_tumors <- categorization_list$mean_similarity_cell_line_to_k_nearest_tumors
  mean_similarity_tumor_to_k_nearest_tumors <- categorization_list$mean_similarity_tumor_to_k_nearest_tumors
  
  num_cell_lines <- length(cell_line_ids)
  num_tumors <- length(tumor_ids)
  
  cell_lines_and_tumors.col <- c(rep(tumor_color, num_cell_lines), rep(cell_line_color, num_tumors))
  cell_lines_and_tumors.pch <- c(rep(tumor_shape, num_cell_lines), rep(cell_line_shape, num_tumors))
  
  if(trim_cell_line_names) {
    cell_line_ids <- sapply(cell_line_ids, return_first_part)
  }  
  
  #plot(isomdsfit$points,col=cell_lines_and_tumors.col, pch=cell_lines_and_tumors.pch, xlab="Coordinate 1", ylab="Coordinate 2" ,main="Weighted Correlation, including low-level CNAs")
  #text(x=isomdsfit$points[1:num_cell_lines,1], y=isomdsfit$points[1:num_cell_lines,2]+0.025, labels=cell_line_ids, cex=0.6)
  
  dataframe_for_ggplot <- as.data.frame(cbind(isomdsfit$points[,1],isomdsfit$points[,2]))
  colnames(dataframe_for_ggplot) <- c("Coordinate1", "Coordinate2")
  dataframe_for_ggplot$Color_for_Sample_Type <- c(rep(cell_line_color,num_cell_lines), rep(tumor_color,num_tumors))
  dataframe_for_ggplot$Shape_for_Sample_Type  <- c(rep(cell_line_shape,num_cell_lines), rep(tumor_shape,num_tumors))
  
  if(use_gradient) {
    tmp <- map_mean_similarity_to_gradient(
      mean_similarity_cell_line_to_k_nearest_tumors=mean_similarity_cell_line_to_k_nearest_tumors,
      mean_similarity_tumor_to_k_nearest_tumors=mean_similarity_tumor_to_k_nearest_tumors,
      col1=cell_line_color,
      col2=tumor_color, 
      numshades=100)
  } else {
    tmp <- rep(cell_line_color, num_cell_lines)
  }
  
  # Use gradients for cell lines and tumor_color for all tumors
  dataframe_for_ggplot$Color_for_Sample_Type_plus_gradient_by_similarity <- c(tmp, rep(tumor_color, num_tumors))

  dataframe_for_ggplot$Labels_for_plot_tumors_blanked  <- c(rep(cell_line_ids),rep("",num_tumors))
  dataframe_for_ggplot$Size_for_Sample_Type  <- c(rep(3, num_cell_lines), rep(1,num_tumors))
  
  ggplot(as.data.frame(dataframe_for_ggplot), aes_string(x="Coordinate1", y="Coordinate2")) + 
    geom_point(colour=dataframe_for_ggplot$Color_for_Sample_Type_plus_gradient_by_similarity, size = dataframe_for_ggplot$Size_for_Sample_Type) + 
    geom_text(label= dataframe_for_ggplot$Labels_for_plot_tumors_blanked) + 
    theme_minimal()
  
  # IGNORE High level only
  #freq_alt_samplewise_cna_high_level_only <- apply(composite_CNA_high_level_only,2,compute_freq_alt)
  #plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
  #text(freq_alt_samplewise_mut[1:length(cell_line_ids)],freq_alt_samplewise_cna_high_level_only[1:length(cell_line_ids)],labels=cell_line_ids,cex=0.6)
  #legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))
  
  # IGNORE High level only
  # isomdsfit2 <- isoMDS(1- cor_weighted_high_level_only,k=2)
  # plot(isomdsfit2$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, excluding low-level CNAs")
  # text(x=isomdsfit2$points[1:num_cell_lines,1], y=isomdsfit2$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)
  # legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))
}
