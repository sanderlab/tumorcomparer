#' Make Balloon Plot for Results 
#' 
#' @param 
#' @param 
#' 
#' @return a ggplot2 object 
#' 
#' @concept tumorcomparer 
#' @export
#' 
#' @importFrom reshape2 melt 
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle xlab ylab labs scale_x_discrete theme_bw theme 
plot_balloon_plot <- function(dat, title) {
  p <- ggplot(dat, aes(x=factor(variable), y=factor(Cell_Line_Name), color=factor(variable), size=value)) + 
    geom_point() + 
    geom_text(aes(label=value, x= as.numeric(factor(variable)) + 0.25), alpha=1.0, size=3) + 
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

