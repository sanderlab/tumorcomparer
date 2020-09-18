#' Make Balloon Plot from MTC Data structure 
#' 
#' @param mtc (NOTE: users should not make dataset; the dataset with results from manuscript 
#'  is included in the package). a data.frame with the columns: 
#'  
#'  * Cell_Line_Name: Name of the cell line 
#'  * Cell_Line_Cancer_Type: Name of the cancer type of the cell line
#'  * Tumor_Cancer_Type: TCGA abbreviation of cancer type for tumor samples
#'  * Mean_Similarity_To_Tumors_AVG: TODO 
#'  * Mean_Similarity_To_Tumors_MUT: TODO
#'  * Mean_Similarity_To_Tumors_CNA: TODO
#'  * Mean_Similarity_To_Tumors_EXP: TODO
#'  * AVGSIM_Zscores; MUTSIM_Zscores: TODO
#'  * CNASIM_Zscores; EXPSIM_Zscores: TODO
#'  * AVGSIM_Percentile_Ranks: TODO
#'  * MUTSIM_Percentile_Ranks: % Rank by Mutation
#'  * CNASIM_Percentile_Ranks: % Rank by Copy Number
#'  * EXPSIM_Percentile_Ranks: % Rank by Expression
#'  * Categorization: TODO  
#'  * AVGSIM_Zscores_wrt_Tumors: TODO
#'  * MUTSIM_Zscores_wrt_Tumors: TODO
#'  * CNASIM_Zscores_wrt_Tumors: TODO 
#'  * EXPSIM_Zscores_wrt_Tumors: TODO 
#'  * Mean_Similarity_To_All_Tumors_AVG: TODO 
#'  * Mean_Similarity_To_All_Tumors_MUT: TODO
#'  * Mean_Similarity_To_All_Tumors_CNA: TODO 
#'  * Mean_Similarity_To_All_Tumors_EXP: TODO
#'  * Average_Of_Percentile_Ranks: TODO
#'  * Rank_of_Average_Of_Percentile_Ranks: % Rank by Avg % Ranks 
#'  
#' @param cancer_type a cancer type abbreviation found in the MTC column Cell_Line_Cancer_Type
#' @param melt_data boolean, whether to apply reshape2::melt() to data.frame
#' 
#' @description Code to convert the processed output of TC over a pan-cancer dataset balloon plots
#'   for a single cancer type ("current_cancer_type") MTC - The subset of TC's output on a
#'   pan-cancer dataest which contains results for matched cancer type, i.e. for each cell line, it
#'   rreports the results when compared to tumors of the same cancer type So, for example, if we run
#'   TC on 10 cancer types on a pan-cancer type, we get 10 sets of results, with each set containing
#'   the results of comparing one of the 10 cancer types to all cell lines From this, MTC extracts
#'   only the results where the cell line and tumors are of the same cancer type, using the other
#'   cell lines to derive ranks We use ranks to compare across cancer types, since the weights and
#'   score distributions are very different across data types and  cancer types
#'   
#' @return a ggplot object
#' 
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle xlab ylab element_blank element_line
#' @importFrom dplyr mutate 
#' @importFrom reshape2 melt 
#' @importFrom stats reorder
#' @importFrom magrittr %>%
#' 
#' @export 
make_balloon_plot_data_from_mtc <- function(mtc, cancer_type, melt_data=TRUE) {
  # NOTE: It is assumed that the mtc data.frame will contain all these columns since they are the results for the publication
  selected_columns <- c("Cell_Line_Name", "MUTSIM_Percentile_Ranks", "CNASIM_Percentile_Ranks", "EXPSIM_Percentile_Ranks", "Rank_of_Average_Of_Percentile_Ranks")
  
  heatmap_mat <- mtc[which(mtc$Cell_Line_Cancer_Type == cancer_type), selected_columns]
  heatmap_mat[,-1] <- apply(heatmap_mat[,-1], 2, as.numeric)
  heatmap_mat[,-1] <- round(heatmap_mat[,-1], digits=2)
  colnames(heatmap_mat)[2:4] <- c("MUTSIM_Ranks", "CNASIM_Ranks", "EXPSIM_Ranks")
  colnames(heatmap_mat)[5] <- "Combined_Score_Ranks"
  
  if(melt_data) {
    df <- melt(heatmap_mat, id.vars = "Cell_Line_Name")
    df <- transform(df, Cell_Line_Name=reorder(Cell_Line_Name, value))     
  } else {
    df <- heatmap_mat
    df <- df[order(-df$Combined_Score_Ranks), ]
  }
  
  return(df)
}
