#' Return first part of ID for CCLE cell line names 
#' 
#' @param x an ID (character string), with fields separated by "_", e.g. "IGROV1_OVARY"
#'
#' @return the part of the string before the first "_", e.g. IGROV1
#' 
#' @details This is only for CCLE (or similar) IDs and is used for plotting
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export 
return_first_part <- function(x) {
  return(strsplit(x,"_")[[1]][1])
}