#' Rescale between 0-1
#' 
#' @description  Map the distance matrices for all data types to 0-1 
#' 
#' @param x a numeric vector of values 
#' 
#' @return a rescaled numeric vector of values  
#' 
#' @export
convert_to_0_to_1_using_xminusmin_by_maxminusmin <- function(x) { (x-min(x))/(max(x)-min(x)) }
