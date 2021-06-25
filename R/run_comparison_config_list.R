run_comparison_config_list <- function(config_list) {
  
  #### checking input data ####
  
  if(length(config_list) < 1) {
    stop("ERROR: At least one data type: mut, cna, exp must be provided for available_data_types")  
  }
  
  
  if(sum(sapply(config_list, function(x) {x$data_type_weight})) != 1) {
    stop("ERROR: Sum of *_data_type_weights must sum up to 1")  
  }
  
  #### runing comparison ####
  
  calculated_datsets_comparison <- lapply(config_list, function(x) {
    
    generate_composite_mat_and_gene_weights(
      default_weight = x$default_weight,
      tumor_file = x$tumor_file,
      cell_line_file = x$cell_line_file,
      known_cancer_gene_weights_file = x$known_cancer_gene_weights_file,
      cancer_specific_gene_weights_file = x$cancer_specific_gene_weights_file)
    
  })
  
  #### harmonization of the datasets ####
  
  combined_samples <- Reduce(intersect, lapply(calculated_datsets_comparison, function(x) {colnames(x$dist_mat)}))
  combined_tumor_ids <- Reduce(intersect, lapply(calculated_datsets_comparison, function(x) {x$tumor_ids}))
  combined_cell_line_ids <- Reduce(intersect, lapply(calculated_datsets_comparison, function(x) {x$cell_line_ids}))
  
  
  # CALCULATE COMBINED_DIST AND ISOMDS ----
  
  combined_dist <- Reduce('+',
         lapply(config_list, function(x) {
           x$data_type_weight * calculated_datsets_comparison[[x$dataset_name]]$dist_mat[combined_samples, combined_samples]
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
    dist_mat_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$dist_mat}),
    composite_mat_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$composite_mat}),
    gene_weights_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$gene_weights}),
    isomdsfit = isomdsfit,
    isomdsfit_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$isomdsfit}),
    cell_line_ids = combined_cell_line_ids,
    tumor_ids = combined_tumor_ids,
    known_cancer_gene_weights_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$known_cancer_gene_weights}),
    cancer_specific_gene_weights_by_data_type = lapply(calculated_datsets_comparison, function(x) {x$cancer_specific_gene_weights})
  )
  
  return(results)
  
}