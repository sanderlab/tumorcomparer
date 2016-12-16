#' 
return_first_part <- function(x) {
  return(strsplit(x,"_")[[1]][1])
}