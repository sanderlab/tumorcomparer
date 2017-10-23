return_fourth_part_of_TCGA_ID <- function(x){
  return(strsplit(as.character(x),"-")[[1]][4])
}

return_first_part_Bar <- function(x){
  return(strsplit(x,"\\|")[[1]][1])
}

return_3_field_TCGA_ID <- function(x){
  return( paste(strsplit((x),"-")[[1]][1],str_split((x),"-")[[1]][2],str_split((x),"-")[[1]][3],sep="-"))
}

return_first_part <- function(x){ return(strsplit(x,"_")[[1]][1]) }

return_second_part <- function(x){ return(strsplit(x,"_")[[1]][2]) }
