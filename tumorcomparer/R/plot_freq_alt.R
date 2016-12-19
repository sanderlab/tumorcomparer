#' Plot mutation and copy number fractions 
#' 
#' @param num_cell_lines number of cell lines from FIX FUNCTION X
#' @param num_tumors number of cell lines from FIX FUNCTION Y
#' @param composite_mut a composite (both tumor and cell line information) matrix with only mutation information
#' @param composite_cna a composite (both tumor and cell line information) matrix with only copy number information
#' @param tumor_color a color for tumor points (DEFAULT: orange)
#' @param cell_line_color a color for tumor points (DEFAULT: blue)
#' @param tumor_shape an integer for an R plot PCH symbol
#' @param cell_line_shape an integer for an R plot PCH symbol
#' 
#' @return Nothing is returned 
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
plot_freq_alt <- function(num_cell_lines, num_tumors, composite_mut, composite_cna, tumor_color="orange", cell_line_color="blue", tumor_shape=17, cell_line_shape=20) {
  #num_cell_lines <- length(cell_lines_with_both_MUT_and_CNA)
  #num_tumors <- length(tumors_with_both_MUT_and_CNA)
  
  cell_lines_and_tumors.col <- c(rep(tumor_color, num_cell_lines), rep(cell_line_color, num_tumors))
  cell_lines_and_tumors.pch <- c(rep(tumor_shape, num_cell_lines),rep(cell_line_shape, num_tumors))
  
  freq_alt_samplewise_mut <- apply(composite_mut, 2, compute_freq_alt)
  freq_alt_samplewise_cna <- apply(composite_cna, 2, compute_freq_alt)
  
  plot(freq_alt_samplewise_mut, freq_alt_samplewise_cna, col=cell_lines_and_tumors.col, pch=cell_lines_and_tumors.pch, xlim=c(0,1), ylim=c(0,1), xlab="Fraction Genes Mutated", ylab="Fraction genes copy number altered",main="Using all CNAs")
  text(freq_alt_samplewise_mut[1:length(cell_line_ids)], freq_alt_samplewise_cna[1:length(cell_line_ids)], labels=cell_line_ids,cex=0.6)
  
  #freq_alt_samplewise_cna_high_level_only <- apply(composite_cna_high_level_only,2,compute_freq_alt)
  #plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
  #text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_lines_with_both_MUT_and_CNA,cex=0.6)
}

