#' Run a comparison between two cohorts (e.g. cell lines and tumors) based on specified data types by configuration list
#' 
#' @param config_list list which sepcifies comparison datasets and their parameters for each data type (mutation, CNV, expression, etc.).
#'    The structure of the config_list is list(mut = list(), cnv = list(), exp = list(), ...). The list of each dataset must contain the following arguments
#'    list(dataset_name, 
#'         data_type_weight, 
#'         default_weight, 
#'         tumor_file, 
#'         cell_line_file, 
#'         known_cancer_gene_weights_file, 
#'         cancer_specific_gene_weights_file)
#'    dataset_name: short name of comparison data type
#'    data_type_weight: a numeric weight for the data type (NOTE: data type weights must sum to 1)
#'    default_weight: default (background) weight for specified data type (EXP) (DEFAULT: 0.01)
#'    tumor_file: path to a dataset file which contain 
#'    cell_line_file:
#'    known_cancer_gene_weights_file:
#'    cancer_specific_gene_weights_file:
#' @param remove_errored_dataset_comparisons will skip the data types which can't be compared for technical reasons when set to TRUE (Default: FALSE)
#' @param gene_list a vector of HGNC gene symbols to run comparison only for the specified genes (Default: NULL)
#' 
#' @return FIXME
#' 
#' @export 
run_comparison_config_list <- function(config_list, gene_list = NULL, remove_errored_dataset_comparisons=FALSE) {
  
  #### checking input data ####
  
  if(length(config_list) < 1) {
    stop("ERROR: At least one data type: mut, cna, exp must be provided for available_data_types")  
  }
  
  
  if(sum(sapply(config_list, function(x) {x$data_type_weight})) != 1) {
    stop("ERROR: Sum of *_data_type_weights must sum up to 1")  
  }
  
  #### running comparison ####
  
  calculated_datasets_comparison <- lapply(config_list, function(x) {
    
    tumor <- read.table(x$tumor_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    
    if(sum(colSums(tumor)) == 0) {
      
      if(remove_errored_dataset_comparisons) {
        cat(paste0("INFO: skipping ", x$dataset_name, " data, zero values only", "\n"))
        
        NULL
        
      } else {
        stop(paste0("ERROR: ", x$dataset_name, " data does not have any values higher than zero"))  
      }
      
    } else {
      
      generate_composite_mat_and_gene_weights(
        default_weight = x$default_weight,
        tumor_file = x$tumor_file,
        cell_line_file = x$cell_line_file,
        known_cancer_gene_weights_file = x$known_cancer_gene_weights_file,
        cancer_specific_gene_weights_file = x$cancer_specific_gene_weights_file,
        gene_list
        )
      
    }
    
  })
  
  ## filtering skipped datasets
  calculated_datasets_comparison <- calculated_datasets_comparison[which(!sapply(calculated_datasets_comparison, is.null))]
  
  config_list <- config_list[names(calculated_datasets_comparison)]
  
  #### harmonization of the datasets ####
  
  combined_samples <- Reduce(intersect, lapply(calculated_datasets_comparison, function(x) {colnames(x$dist_mat)}))
  combined_tumor_ids <- Reduce(intersect, lapply(calculated_datasets_comparison, function(x) {x$tumor_ids}))
  combined_cell_line_ids <- Reduce(intersect, lapply(calculated_datasets_comparison, function(x) {x$cell_line_ids}))
  
  
  # CALCULATE COMBINED_DIST AND ISOMDS ----
  
  combined_dist <- Reduce('+',
         lapply(config_list, function(x) {
           x$data_type_weight * calculated_datasets_comparison[[x$dataset_name]]$dist_mat[combined_samples, combined_samples]
         })
         )
  
  # RUN ISOMDS ----
  isomdsfit <-  isoMDS(combined_dist, k=2)
  
  # CHECK ----
  if(length(combined_cell_line_ids) == 0) {
    stop("ERROR: Result validation error. Please check that tumor and cell line input files are set correctly.")  
  }
  
  if(length(combined_tumor_ids) == 0) {
    stop("ERROR: Result validation error. Please check that tumor and cell line input files are set correctly.")  
  }
  
  # MERGE RESULTS ----
  results <- list(
    dist_mat = combined_dist,
    dist_mat_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$dist_mat}),
    composite_mat_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$composite_mat}),
    gene_weights_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$gene_weights}),
    isomdsfit = isomdsfit,
    isomdsfit_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$isomdsfit}),
    cell_line_ids = combined_cell_line_ids,
    tumor_ids = combined_tumor_ids,
    known_cancer_gene_weights_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$known_cancer_gene_weights}),
    cancer_specific_gene_weights_by_data_type = lapply(calculated_datasets_comparison, function(x) {x$cancer_specific_gene_weights}),
    calculated_data_types = names(calculated_datasets_comparison)
  )
  
  return(results)
  
}