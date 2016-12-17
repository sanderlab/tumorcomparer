#'Return first part of ID/string
#' Input: @x an ID (character string), with fields separated by "_", e.g. "IGROV1_OVARY"
#' Output: The part of the string before the first "_", e.g. IGROV1
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (augustin@mail.nih.gov)
#'
#' 
#'
#' @concept tumorcomparer
#' @export 
return_first_part <- function(x) {
  return(strsplit(x,"_")[[1]][1])
}