#' Plot mutation and copy number fractions 
#' 
#' @param comparison_list the results of run_comparison() (See: run_comparison())
#' @param tumor_color a color for tumor points (DEFAULT: blue)
#' @param cell_line_color a color for tumor points (DEFAULT: orange)
#' @param tumor_shape an integer for an R plot PCH symbol
#' @param cell_line_shape an integer for an R plot PCH symbol
#' @param trim_cell_line_names a boolean whether to trim the the cell lines; 
#'   this is optional and used for long cell line names in CCLE format (i.e. CELLLINE_TISSUE)
#' @param max_xlim a value between 0 and 1 the maximum value to for x-axis
#' 
#' @return Nothing is returned 
#' 
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom ggplot2 ggplot aes geom_point geom_text theme_minimal
plot_freq_alt <- function(comparison_list, tumor_color="blue", cell_line_color="orange", tumor_shape=17, cell_line_shape=20, trim_cell_line_names=TRUE, max_xlim=1) {

  cell_lines_with_both_MUT_and_CNA <- comparison_list$cell_lines_with_both_MUT_and_CNA
  tumors_with_both_MUT_and_CNA <- comparison_list$tumors_with_both_MUT_and_CNA
  composite_MUT <- comparison_list$composite_MUT
  composite_CNA <- comparison_list$composite_CNA
  isomdsfit <- comparison_list$isomdsfit
  
  num_cell_lines <- length(cell_lines_with_both_MUT_and_CNA)
  num_tumors <- length(tumors_with_both_MUT_and_CNA)
  
  cell_lines_and_tumors.col <- c(rep(tumor_color, num_cell_lines), rep(cell_line_color, num_tumors))
  cell_lines_and_tumors.pch <- c(rep(tumor_shape, num_cell_lines),rep(cell_line_shape, num_tumors))
  
  freq_alt_samplewise_mut <- apply(composite_MUT, 2, compute_freq_alt)
  freq_alt_samplewise_cna <- apply(composite_CNA, 2, compute_freq_alt)
  
  if(trim_cell_line_names) {
    cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, return_first_part)
  } else {
    cell_line_ids <- cell_lines_with_both_MUT_and_CNA
  }
  
  # Build necessary data.frame
  dataframe_for_ggplot <- data.frame(
    Coordinate1=isomdsfit$points[,1], 
    Coordinate2=isomdsfit$points[,2], 
    Color_for_Sample_Type=c(rep(cell_line_color,num_cell_lines), rep(tumor_color,num_tumors)),
    Shape_for_Sample_Type=c(rep(cell_line_shape,num_cell_lines), rep(tumor_shape,num_tumors)),
    Fraction_Genes_Mutated=c(freq_alt_samplewise_mut),
    Fraction_Genes_Copy_Number_Altered=c(freq_alt_samplewise_cna),
    Labels_for_plot_tumors_blanked=c(rep(cell_line_ids),rep("",num_tumors)),
    Size_for_Sample_Type=c(rep(3, num_cell_lines), rep(1,num_tumors)))
  
  ggplot(dataframe_for_ggplot, aes_string(x="Fraction_Genes_Mutated", y="Fraction_Genes_Copy_Number_Altered")) +
    geom_point(colour= dataframe_for_ggplot$Color_for_Sample_Type, size = dataframe_for_ggplot$Size_for_Sample_Type) + 
    geom_text(label= dataframe_for_ggplot$Labels_for_plot_tumors_blanked) + 
    theme_minimal()
  
  # plot(freq_alt_samplewise_mut, freq_alt_samplewise_cna, col=cell_lines_and_tumors.col, pch=cell_lines_and_tumors.pch, xlim=c(0,max_xlim), ylim=c(0,1), xlab="Fraction Genes Mutated", ylab="Fraction genes copy number altered",main="Frequency of Alterations in Cell Lines and Tumors")
  # text(freq_alt_samplewise_mut[1:length(cell_line_ids)], freq_alt_samplewise_cna[1:length(cell_line_ids)], labels=cell_line_ids,cex=0.6)
  
  #freq_alt_samplewise_cna_high_level_only <- apply(composite_CNA_high_level_only,2,compute_freq_alt)
  #plot(freq_alt_samplewise_mut,freq_alt_samplewise_cna_high_level_only,col=cell_lines_and_tumors.col,xlab="Fraction Genes Mutated",ylab="Fraction genes copy number altered",pch=cell_lines_and_tumors.pch,ylim=c(0,1),xlim=c(0,1),main="Using high-level CNAs only")
  #text(freq_alt_samplewise_mut[1:length(cell_lines_with_both_MUT_and_CNA)],freq_alt_samplewise_cna_high_level_only[1:length(cell_lines_with_both_MUT_and_CNA)],labels=cell_lines_with_both_MUT_and_CNA,cex=0.6)
}

