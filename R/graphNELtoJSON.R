#' Convert R graphNEL object to cytoscape.js JSON.
#' 
#' @param g a graphNEL
#'
#' @concept tumorcomparer
#' @export
#' 
#' @import graph
#' @import jsonlite
graphNELtoJSON <- function (g) 
{
  if (length(nodes(g)) == 0) {
    return("{}")
  }
  vector.count <- 10 * (length(edgeNames(g)) + length(nodes(g)))
  vec <- vector(mode = "character", length = vector.count)
  i <- 1
  vec[i] <- "{\"elements\": {\"nodes\": ["
  i <- i + 1
  nodes <- nodes(g)
  edgeNames <- edgeNames(g)
  edges <- strsplit(edgeNames, "~")
  edgeNames <- sub("~", "->", edgeNames)
  names(edges) <- edgeNames
  noa.names <- names(graph::nodeDataDefaults(g))
  eda.names <- names(graph::edgeDataDefaults(g))
  nodeCount <- length(nodes)
  edgeCount <- length(edgeNames)
  for (n in 1:nodeCount) {
    node <- nodes[n]
    vec[i] <- "{\"data\": "
    i <- i + 1
    nodeList <- list(id = node)
    this.nodes.data <- graph::nodeData(g, node)[[1]]
    if (length(this.nodes.data) > 0) {
      nodeList <- c(nodeList, this.nodes.data)
    }
    nodeList.json <- toJSON(nodeList, auto_unbox = TRUE)
    vec[i] <- nodeList.json
    i <- i + 1
    if (all(c("xPos", "yPos") %in% names(graph::nodeDataDefaults(g)))) {
      position.markup <- sprintf(", \"position\": {\"x\": %f, \"y\": %f}", 
                                 graph::nodeData(g, node, "xPos")[[1]], graph::nodeData(g, 
                                                                                        node, "yPos")[[1]])
      vec[i] <- position.markup
      i <- i + 1
    }
    if (all(c("x", "y") %in% names(graph::nodeDataDefaults(g)))) {
      position.markup <- sprintf(", \"position\": {\"x\": %f, \"y\": %f}", 
                                 graph::nodeData(g, node, "x")[[1]], graph::nodeData(g, 
                                                                                     node, "y")[[1]])
      vec[i] <- position.markup
      i <- i + 1
    }
    if (n != nodeCount) {
      vec[i] <- "},"
      i <- i + 1
    }
  }
  vec[i] <- "}]"
  i <- i + 1
  if (edgeCount > 0) {
    vec[i] <- ", \"edges\": ["
    i <- i + 1
    for (e in seq_len(edgeCount)) {
      vec[i] <- "{\"data\": "
      i <- i + 1
      edgeName <- edgeNames[e]
      edge <- edges[[e]]
      sourceNode <- edge[[1]]
      targetNode <- edge[[2]]
      edgeList <- list(id = edgeName, source = sourceNode, 
                       target = targetNode)
      this.edges.data <- graph::edgeData(g, sourceNode, 
                                         targetNode)[[1]]
      if (length(this.edges.data) > 0) {
        edgeList <- c(edgeList, this.edges.data)
      }
      edgeList.json <- toJSON(edgeList, auto_unbox = TRUE)
      vec[i] <- edgeList.json
      i <- i + 1
      if (e != edgeCount) {
        vec[i] <- "},"
        i <- i + 1
      }
    }
    vec[i] <- "}]"
    i <- i + 1
  }
  vec[i] <- "}"
  i <- i + 1
  vec[i] <- "}"
  vec.trimmed <- vec[which(vec != "")]
  paste0(vec.trimmed, collapse = " ")
}