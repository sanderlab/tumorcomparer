#' Make Balloon Plot for Results 
#' 
#' @param dat data generated from either make_
#' @param title a string to be the title for the plot
#' 
#' @return a ggplot2 object 
#' 
#' @concept tumorcomparer 
#' @export
#' 
#' @seealso [make_balloon_plot_data_from_comparison_result()], 
#' [make_balloon_plot_data_from_mtc()]
#' 
#' @examples 
#' mtc_file <- system.file("extdata", "mtc_results_20200331", "mtc_results_20200331_no_factors.rds", 
#'   package="tumorcomparer")
#' mtc <- readRDS(mtc_file)
#' dat <- make_balloon_plot_data_from_mtc(mtc, "BLCA")
#' plot_balloon_plot(dat, "Title")
#' 
#' @importFrom reshape2 melt 
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle xlab ylab labs scale_x_discrete theme_bw theme 
plot_balloon_plot <- function(dat, title) {
  dat$Cell_Line_Name <- factor(dat$Cell_Line_Name)
  dat$variable <- factor(dat$variable)
  dat$numeric_variable <- as.numeric(dat$variable) + 0.25
  
  p <- ggplot(dat, aes_string(x='variable', y='Cell_Line_Name', color='variable', size='value')) + 
    geom_point() + 
    geom_text(aes_string(label='value', x='numeric_variable'), alpha=1.0, size=3) + 
    ggtitle(title) +
    xlab("Weighted Similarity Ranks By Data Type") + 
    ylab ("Cell Lines") + 
    labs(colour="Color", size="Size") +
    #scale_x_discrete(labels = c('Mutation', 'CNA', 'Expression', 'Average')) +
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

