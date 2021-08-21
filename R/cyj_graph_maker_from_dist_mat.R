#' Creates cytoscape.js JSON graph from distance matrix for defined threshold
#' 
#' @param dist_mat distance matrix calculated by run_comparison_config_list() function
#' @param min_weight numeric value between 0 and 1
#' 
#' @return returns cytoscape.js JSON graph
#'   
#'
#' @author Hakobyan Siras (sirashakobyan@gmail.com)
#'
#' @concept tumorcomparer
#' @export
#'
cyj_graph_maker_from_dist_mat <- function(dist_mat, min_weight) {
  
  
  dist_mat_melted <- reshape2::melt(`is.na<-`(dist_mat, upper.tri(dist_mat, diag = T)), na.rm = TRUE)
  
  g <- graph::ftM2graphNEL(as.matrix(dist_mat_melted[which(dist_mat_melted$value > min_weight),1:2]), edgemode = "directed")
  
  graph::nodeDataDefaults(g, attr = "label") <- "NA"
  graph::nodeDataDefaults(g, attr = "color") <- "NA"
  graph::nodeDataDefaults(g, attr = "xPos") <- 0
  graph::nodeDataDefaults(g, attr = "yPos") <- 0
  
  graph::nodeData(g, attr = "label") <- graph::nodes(g)
  graph::nodeData(g, graph::nodes(g)[grepl("TCGA", graph::nodes(g))], attr = "color") <- "#99ccff"
  graph::nodeData(g, graph::nodes(g)[!grepl("TCGA", graph::nodes(g))], attr = "color") <- "#49d849"
  graph::nodeData(g, attr = "xPos") <- igraph::layout_nicely(igraph::graph_from_graphnel(g))[,1]*100
  graph::nodeData(g, attr = "yPos") <- igraph::layout_nicely(igraph::graph_from_graphnel(g))[,2]*100
  
  graph::edgeDataDefaults(g, attr = "edgeType") <- "pp"
  graph::edgeDataDefaults(g, attr = "dist") <- 0
  
  graph::edgeData(g, from = as.character(dist_mat_melted[which(dist_mat_melted$value > min_weight),1]), 
                  to = as.character(dist_mat_melted[which(dist_mat_melted$value > min_weight),2]), 
                  attr = "dist") <- round(dist_mat_melted[which(dist_mat_melted$value > min_weight),3], digits = 2)
  
  return(graphNELtoJSON(g))
  
}