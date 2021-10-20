#' Make Joy Plots for Distributions for Similarity Distributions
#' 
#' @param comparison_result the results of run_comparison() (See: run_comparison())
#' 
#' @note TODO Abstract code to many datasets; use list (https://wilkelab.org/cowplot/reference/plot_grid.html)
#'   
#' @return a ggplot object
#' 
#' @examples 
#' # Generated using: tumorcomparer::run_comparison() 
#' comparison_result <- readRDS(system.file("test_output", "ov_comparison_result.rds", 
#'   package="tumorcomparer"))
#'   
#' plot_joyplot(comparison_result = comparison_result)
#' 
#' @importFrom ggplot2 ggplot geom_point geom_text ggtitle xlab ylab element_blank element_line margin
#' @importFrom ggridges geom_density_ridges
#' @importFrom cowplot plot_grid ggdraw draw_label
#' @importFrom reshape2 melt 
#' @importFrom stats reorder
#' @importFrom magrittr %>%
#' 
#' @export 
plot_joyplot <- function(comparison_result) {
  # Distance matrix from comparison_result
  dist_mat <- comparison_result$dist_mat
  
  # Distance matrix by data type from comparison_result 
  dist_mat_by_data_type <- comparison_result$dist_mat_by_data_type
  
  tumor_ids <- comparison_result$tumor_ids
  cell_line_ids <- comparison_result$cell_line_ids
  
  mut_dat <- (1-dist_mat_by_data_type$mut[cell_line_ids[order(rowMeans(dist_mat[cell_line_ids, tumor_ids]))],tumor_ids])
  cna_dat <- (1-dist_mat_by_data_type$cna[cell_line_ids[order(rowMeans(dist_mat[cell_line_ids,tumor_ids]))],tumor_ids])
  exp_dat <- (1-dist_mat_by_data_type$exp[cell_line_ids[order(rowMeans(dist_mat[cell_line_ids,tumor_ids]))],tumor_ids])
  
  normalized_dat <- 
    convert_to_0_to_1_using_xminusmin_by_maxminusmin(mut_dat) + 
    convert_to_0_to_1_using_xminusmin_by_maxminusmin(cna_dat) + 
    convert_to_0_to_1_using_xminusmin_by_maxminusmin(exp_dat)
  
  data_type_count <- length(comparison_result$dist_mat_by_data_type)
  normalized_and_averaged <- normalized_dat/data_type_count
  
  # PARAMETERS ----
  xlab_title <- "Weighted Similarity to Tumors"
  ylab_title <- "Cell Lines"
  overall_title <- "Weighted Similarity (Averaged and By Specific Profiling)"
  
  # PLOTS ----
  dat <- melt(cna_dat)
  dat$Var1 <- as.factor(dat$Var1)
  
  joyplot_cna <- ggplot(dat, aes_string(y="Var1", x="value")) + 
    geom_density_ridges(alpha=0.5) + 
    labs(title = "By Copy Number Aberrations") + 
    xlab(xlab_title) + 
    ylab(ylab_title) +
    theme_bw()
  joyplot_cna
  
  dat <- melt(mut_dat)
  dat$Var1 <- as.factor(dat$Var1)
  joyplot_mut <- ggplot(dat, aes_string(y="Var1", x="value")) + 
    geom_density_ridges(alpha=0.5) + 
    xlab(xlab_title) + 
    ylab(ylab_title) + 
    labs(title = "By Mutations") +
    theme_bw()
  joyplot_mut
  
  dat <- melt(exp_dat)
  dat$Var1 <- as.factor(dat$Var1)
  joyplot_exp <- ggplot(dat, aes_string(y="Var1", x="value")) + 
    geom_density_ridges(alpha=0.5) + 
    xlab(xlab_title) + 
    ylab(ylab_title) + 
    labs(title = "By Gene Expression") +
    theme_bw()
  joyplot_exp
  
  dat <- melt(normalized_and_averaged[rev(order(rowMeans(normalized_and_averaged))),])
  dat$Var1 <- as.factor(dat$Var1)
  joyplot_avg_after_normalization <- ggplot(dat, aes_string(y="Var1", x="value")) + 
    geom_density_ridges(alpha=0.5) + 
    xlab(xlab_title) + 
    ylab(ylab_title) + 
    labs(title = "Average Weighted Similarity") +
    theme_bw()
  joyplot_avg_after_normalization
  
  #plot_grid(joyplot_avg_after_normalization, joyplot_exp, joyplot_mut, joyplot_cna, ncol=2, labels=c('A', 'B', 'C', 'D'))
  joyplot_composite <- plot_grid(joyplot_avg_after_normalization, joyplot_exp, joyplot_mut, joyplot_cna, ncol=2, labels=c('A', 'B', 'C', 'D'))
  joyplot_composite
  
  ## Overall title 
  title <- ggdraw() + 
    draw_label(
      overall_title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  grid_with_title <- plot_grid(
    title, 
    joyplot_composite,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  grid_with_title
  
  #save_plot("joyplots.png", grid_with_title, base_height = 10, base_width = 10, units = "in")
}