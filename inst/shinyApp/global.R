options(shiny.maxRequestSize = 24*1024^2) # Line exists in server, ui, global

# DEBUG ----
cat("MAX UPLOAD SIZE: ", getOption("shiny.maxRequestSize"), "\n")

# GENERAL ---- 
plot_title_prefix <-  "Mean Similarity to Tumors" 

plotlyModeBarButtonsToRemove <- c(
  "select2d", "sendDataToCloud", "pan2d", "resetScale2d", "hoverClosestCartesian", 
  "hoverCompareCartesian", "lasso2d", "zoomIn2d", "zoomOut2d", "toggleSpikelines"
)

# LOAD DATA ----
mtc_file <- system.file('extdata/mtc_results_20200331/mtc_results_20200331.rds', package="tumorcomparer")
#mtc_file <- system.file('extdata/mtc_results_20200331/mtc_results_20200331_no_factors.rds', package="tumorcomparer")
mtc_dataset <- readRDS(mtc_file)
precomputed_comparisons <- readRDS(system.file('extdata/precomputed_comparisons.rds', package="tumorcomparer"))
selected_geneset_comparions <- readRDS(system.file('extdata/selected_geneset_comparions.rds', package="tumorcomparer"))

mtc_dataset$Cell_Line_Name <- as.character(mtc_dataset$Cell_Line_Name)
#mtc_dataset$Cell_Line_Cancer_Type <- as.character(mtc_dataset$Cell_Line_Cancer_Type)
#mtc_dataset$Tumor_Cancer_Type <- as.character(mtc_dataset$Tumor_Cancer_Type)
mtc_dataset$Mean_Similarity_To_Tumors_AVG <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_AVG))[mtc_dataset$Mean_Similarity_To_Tumors_AVG]
mtc_dataset$Mean_Similarity_To_Tumors_MUT <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_MUT))[mtc_dataset$Mean_Similarity_To_Tumors_MUT]
mtc_dataset$Mean_Similarity_To_Tumors_CNA <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_CNA))[mtc_dataset$Mean_Similarity_To_Tumors_CNA] 
mtc_dataset$Mean_Similarity_To_Tumors_EXP <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_EXP))[mtc_dataset$Mean_Similarity_To_Tumors_EXP]
mtc_dataset$AVGSIM_Zscores <- as.numeric(levels(mtc_dataset$AVGSIM_Zscores))[mtc_dataset$AVGSIM_Zscores]
mtc_dataset$MUTSIM_Zscores <- as.numeric(levels(mtc_dataset$MUTSIM_Zscores))[mtc_dataset$MUTSIM_Zscores]
mtc_dataset$CNASIM_Zscores <- as.numeric(levels(mtc_dataset$CNASIM_Zscores))[mtc_dataset$CNASIM_Zscores]
mtc_dataset$EXPSIM_Zscores <- as.numeric(levels(mtc_dataset$EXPSIM_Zscores))[mtc_dataset$EXPSIM_Zscores]
mtc_dataset$AVGSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$AVGSIM_Percentile_Ranks))[mtc_dataset$AVGSIM_Percentile_Ranks]
mtc_dataset$MUTSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$MUTSIM_Percentile_Ranks))[mtc_dataset$MUTSIM_Percentile_Ranks]
mtc_dataset$CNASIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$CNASIM_Percentile_Ranks))[mtc_dataset$CNASIM_Percentile_Ranks]
mtc_dataset$EXPSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$EXPSIM_Percentile_Ranks))[mtc_dataset$EXPSIM_Percentile_Ranks]
#FIXME: mtc_dataset$Categorization <- as.character(mtc_dataset$Categorization)
mtc_dataset$AVGSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$AVGSIM_Zscores_wrt_Tumors))[mtc_dataset$AVGSIM_Zscores_wrt_Tumors]
mtc_dataset$MUTSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$MUTSIM_Zscores_wrt_Tumors))[mtc_dataset$MUTSIM_Zscores_wrt_Tumors]
mtc_dataset$CNASIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$CNASIM_Zscores_wrt_Tumors))[mtc_dataset$CNASIM_Zscores_wrt_Tumors]
mtc_dataset$EXPSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$EXPSIM_Zscores_wrt_Tumors))[mtc_dataset$EXPSIM_Zscores_wrt_Tumors]
mtc_dataset$Mean_Similarity_To_All_Tumors_AVG <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_AVG))[mtc_dataset$Mean_Similarity_To_All_Tumors_AVG]
mtc_dataset$Mean_Similarity_To_All_Tumors_MUT <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_MUT))[mtc_dataset$Mean_Similarity_To_All_Tumors_MUT]
mtc_dataset$Mean_Similarity_To_All_Tumors_CNA <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_CNA))[mtc_dataset$Mean_Similarity_To_All_Tumors_CNA]
mtc_dataset$Mean_Similarity_To_All_Tumors_EXP <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_EXP))[mtc_dataset$Mean_Similarity_To_All_Tumors_EXP]
mtc_dataset$Average_Of_Percentile_Ranks <- as.numeric(levels(mtc_dataset$Average_Of_Percentile_Ranks))[mtc_dataset$Average_Of_Percentile_Ranks]
mtc_dataset$Rank_of_Average_Of_Percentile_Ranks <- as.numeric(levels(mtc_dataset$Rank_of_Average_Of_Percentile_Ranks))[mtc_dataset$Rank_of_Average_Of_Percentile_Ranks]

#saveRDS(mtc_dataset, 'inst/extdata/mtc_results_20200331/mtc_results_20200331_no_factors.rds')

# Column name mapping
mtc_selected_columns <- c(
  "Cell Line"="Cell_Line_Name", 
  #"Cell Line Type"="Cell_Line_Cancer_Type", 
  #"Tumor Type"="Tumor_Cancer_Type",
  "% Rank by Mutation"="MUTSIM_Percentile_Ranks", 
  "% Rank by Copy Number"="CNASIM_Percentile_Ranks", 
  "% Rank by Expression"="EXPSIM_Percentile_Ranks",
  "% Rank by Avg % Ranks"="Rank_of_Average_Of_Percentile_Ranks")

comparison_result_columns <- c(
  "Cell Line"="Cell_Line_Name", 
  #"Cell Line Type"="Cell_Line_Cancer_Type", 
  #"Tumor Type"="Tumor_Cancer_Type",
  "Mutation Score"="mut_score", 
  "Copy Number Score"="cna_score", 
  "Expression Score"="exp_score",
  "Combined Score"="combined_score") 

tcgaTypes <- c(
  "Adrenocortical Carcinoma"="ACC",
  "Bladder Urothelial Carcinoma"="BLCA",
  "Breast Invasive Carcinoma"="BRCA",
  "Cervical Squamous Cell Carcinoma"="CESC",
  "Cholangiocarcinoma"="CHOL",
  "Colon Adenocarcinoma"="COAD",
  "Diffuse Large B-Cell Lymphoma"="DLBC",
  "Esophageal Adenocarcinoma"="ESCA",
  "Glioblastoma Multiforme"="GBM", 
  "Head and Neck Squamous Cell Carcinoma"="HNSC",
  "Kidney Renal Clear Cell Carcinoma"="KIRC",
  "Acute Myeloid Leukemia"="LAML",
  "Low-Grade Glioma"="LGG",
  "Liver Hepatocellular Carcinoma"="LIHC",
  "Lung adenocarcinoma"="LUAD",   
  "Lung squamous cell"="LUSC", 
  "Mesothelioma"="MESO",
  "Ovarian"="OV",
  "Pancreatic Adenocarcinoma"="PAAD",
  "Prostate Adenocarcinoma"="PRAD",
  "Rectal Adenocarcinoma"="READ",
  "Cutaneous Melanoma"="SKCM",
  "Stomach Adenocarcinoma"="STAD",
  "Thyroid Cancer"="THCA",
  "Endometrial Carcinoma"="UCEC"
)
tcgaTypes <- tcgaTypes[tcgaTypes %in% as.character(unique(mtc_dataset$Tumor_Cancer_Type))]


genesets <- c("Most Variable Genes", "Cell Cycle", "HIPPO", "MYC", "NOTCH", "PI3K", "RTK RAS", "TGF-Beta", "WNT")


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

cyj_graph_maker_from_dist_mat <- function(dist_mat, min_weight) {
  
  
  dist_mat_melted <- reshape2::melt(`is.na<-`(dist_mat, upper.tri(dist_mat, diag = T)), na.rm = TRUE)
  
  g <- graph::ftM2graphNEL(as.matrix(dist_mat_melted[which(dist_mat_melted$value > min_weight),1:2]), edgemode = "directed")
  
  graph::nodeDataDefaults(g, attr = "label") <- "NA"
  graph::nodeDataDefaults(g, attr = "color") <- "NA"
  
  graph::nodeData(g, attr = "label") <- graph::nodes(g)
  graph::nodeData(g, graph::nodes(g)[grepl("TCGA", graph::nodes(g))], attr = "color") <- "#99ccff"
  graph::nodeData(g, graph::nodes(g)[!grepl("TCGA", graph::nodes(g))], attr = "color") <- "#49d849"
  
  graph::edgeDataDefaults(g, attr = "edgeType") <- "pp"
  graph::edgeDataDefaults(g, attr = "dist") <- 0
  
  graph::edgeData(g, from = as.character(dist_mat_melted[which(dist_mat_melted$value > min_weight),1]), 
                  to = as.character(dist_mat_melted[which(dist_mat_melted$value > min_weight),2]), 
                  attr = "dist") <- round(dist_mat_melted[which(dist_mat_melted$value > min_weight),3], digits = 2)
  
  return(cyjShiny::graphNELtoJSON(g))
  
}