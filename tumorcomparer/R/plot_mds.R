#' Run a comparison between 
#' 
#' @param cell_lines_with_both_MUT_and_CNA a vector of cell line IDs/names with 
#'   both mutation (MUT) and copy number alteration (CNA) information (see run_comparison)
#' @param tumors_with_both_MUT_and_CNA a vector of tumor IDs with both mutation 
#'   (MUT) and copy number alteration (CNA) information
#' @param tumor_color a color for tumor points (DEFAULT: orange)
#' @param cell_line_color a color for tumor points (DEFAULT: blue)
#' @param tumor_shape an integer for an R plot PCH symbol (DEFAULT: 17)
#' @param cell_line_shape an integer for an R plot PCH symbol (DEFAULT: 20)
#' @param cell_line_ids a vector of strings with cell line names. This will be 
#'   taken from cell_lines_with_both_MUT_and_CNA and assumed to be in the CCLE format 
#'   (i.e. CELLLINE_TISSUE); (DEFAULT: NULL, IDs will be taken from parameter 
#'   cell_lines_with_both_MUT_and_CNA)
#' @param dist_mat a matrix of distances 
#' 
#' @return Nothing is returned 
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom graphics plot text
plot_mds <- function(cell_lines_with_both_MUT_and_CNA, 
                     tumors_with_both_MUT_and_CNA, 
                     dist_mat, 
                     cell_line_ids=NULL,
                     tumor_color="orange", 
                     cell_line_color="blue", 
                     tumor_shape=17, 
                     cell_line_shape=20) {
  # NOTE: Dist should be 
  #dist_mat <- 1-cor_weighted
  
  num_cell_lines <- length(cell_lines_with_both_MUT_and_CNA)
  num_tumors <- length(tumors_with_both_MUT_and_CNA)
  
  cell_lines_and_tumors.col <- c(rep(tumor_color, num_cell_lines), rep(cell_line_color, num_tumors))
  cell_lines_and_tumors.pch <- c(rep(tumor_shape, num_cell_lines), rep(cell_line_shape, num_tumors))
  
  if(!is.null(cell_line_ids)) {
    cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, return_first_part)
  }
  
  isomdsfit <- isoMDS(dist_mat, k=2)
  plot(isomdsfit$points,col=cell_lines_and_tumors.col, pch=cell_lines_and_tumors.pch, xlab="Coordinate 1", ylab="Coordinate 2" ,main="Weighted Correlation, including low-level CNAs")
  text(x=isomdsfit$points[1:num_cell_lines,1], y=isomdsfit$points[1:num_cell_lines,2]+0.025, labels=cell_line_ids, cex=0.6)
  
  # IGNORE High level only
  #freq_alt_samplewise_cna_high_level_only <- apply(composite_CNA_high_level_only,2,compute_freq_alt)
  #plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
  #text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_line_ids,cex=0.6)
  #legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))
  
  
  # IGNORE High level only
  # isomdsfit2 <- isoMDS(1- cor_weighted_high_level_only,k=2)
  # plot(isomdsfit2$points,col=cell_lines_and_tumors.col,pch=cell_lines_and_tumors.pch,xlab="Coordinate 1", ylab="Coordinate 2",main="Weighted Correlation, excluding low-level CNAs")
  # text(x=isomdsfit2$points[1:num_cell_lines,1], y=isomdsfit2$points[1:num_cell_lines,2]+0.025,labels=cell_line_ids,cex=0.6)
  # legend("bottomright",c("Tumors", "Cell Lines"), pch=c(20, 17), cex=.8, col=c("blue", "orange"))
}