#' Transforms data frame of balloon plot into data frame with weights for each data types in each column
#' 
#' @param plot_data output file of make_balloon_plot_data_from_comparison_result() function
#' 
#' @return returns a data frame with weights for each data types in each column
#'   
#'
#' @author Hakobyan Siras (sirashakobyan@gmail.com)
#' @export
#'
#' @concept tumorcomparer
#' @importFrom stats setNames
ballon_plot_data_to_result_table <- function(plot_data) {
  
  comp_table_colnames <- setNames(nm = c("mut_score", "cna_score", "exp_score", "combined_score"), 
                                  object = c("% Rank by Mutation", "% Rank by Copy Number", "% Rank by Expression", "% Rank by Avg % Ranks"))
  
  plot_data <- plot_data$plot_data
  
  plot_data$Cell_Line_Name <- as.character(plot_data$Cell_Line_Name)
  plot_data$variable <- as.character(plot_data$variable)
  
  splitted_plot_data <- split(x = plot_data[,-c(1, 2)], f = plot_data$variable)
  
  merged_df <- as.data.frame(Reduce(cbind, splitted_plot_data))
  
  merged_df <- cbind(split(x = plot_data[,-c(2,3)], f = plot_data$variable)[[1]], merged_df)
  
  colnames(merged_df) <- c("Cell Line", comp_table_colnames[names(splitted_plot_data)] )
  
  merged_df <- merged_df[,c("Cell Line", comp_table_colnames[which(comp_table_colnames %in% colnames(merged_df)[-1])])]
  
  return(merged_df)
}