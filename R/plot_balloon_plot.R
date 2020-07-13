#' Make Balloon Plot for Results 
#' 
#' @param mtc a data.frame with results; check inst/extdata/mtc_results_20200331.rds for format 
#' @param cancer_type a TCGA acronym for a cancer type (e.g., BLCA for bladder cancer)
#' 
#' @return a ggplot2 object 
#' 
#' @note FIXME EXAMPLE #load("data/Results_Matching_Cancer_Types_Results.rda"); #current_cancer_type <- "BLCA"
#' 
#' @concept tumorcomparer 
#' @export
#' 
#' @importFrom reshape2 melt 
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle xlab ylab labs scale_x_discrete theme_bw theme 
plot_balloon_plot <- function(mtc, cancer_type) {
  # TRANSFORMATION ----
  mat_for_heatmap <- mtc[which(mtc$Cell_Line_Cancer_Type == cancer_type),c(1,13:15)]
  mat_for_heatmap[,5] <- as.numeric(mtc[which(mtc$Cell_Line_Cancer_Type == cancer_type),26])
  mat_for_heatmap[,-1] <- apply(mat_for_heatmap[,-1],2,as.numeric)
  mat_for_heatmap[,-1] <- round(mat_for_heatmap[,-1],digits=2)
  mat_for_heatmap[,6] <- round(rowMeans(mat_for_heatmap[,2:4]),digits = 2)
  colnames(mat_for_heatmap)[2:4] <- c("MUTSIM_Ranks","CNASIM_Ranks","EXPSIM_Ranks")
  colnames(mat_for_heatmap)[6] <- "AVGSIM_Scores"
  colnames(mat_for_heatmap)[5] <- "Combined_Score_Ranks"
  mat_for_heatmap <- mat_for_heatmap[,-5]
  df <- reshape2::melt(mat_for_heatmap,id.vars = "Cell_Line_Name")
  
  # PLOT ----
  tmp_df <- df %>% mutate(Cell_Line_Name = reorder(Cell_Line_Name, value, fun=mean)) 
  title <- paste0(cancer_type, " Similarity Ranks by Data Type")
  p <- ggplot(tmp_df, aes(x=factor(variable), y=factor(Cell_Line_Name), color=factor(variable), size=value)) + 
    geom_point() + 
    geom_text(aes(label=value, x= as.numeric(factor(variable)) + 0.25), alpha=1.0, size=3) + 
    ggtitle(title) +
    xlab("") + 
    ylab ("Cell Lines") + 
    labs(colour="Color", size="Size") +
    scale_x_discrete(labels = c('Mutation', 'CNA', 'Expression', 'Average')) +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="none")
  return(p)
  
  #filename <- file.path("balloon_plots", paste0(cancer_type, "_balloon_plot.png"))
  #ggsave(filename, height = 5, width = 6, units = "in") # 7 is probably the minimum
}

# IGNORE ---- 
# Edit legend title and labels
#p + scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))

