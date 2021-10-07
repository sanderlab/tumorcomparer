#' Run a comparison between two cohorts (e.g. cell lines and tumors) based on specified data types by configuration list
#' 
#' @param config_list list which sepcifies comparison datasets and their parameters for each data type (mutation, CNV, expression, etc.).
#'    The structure of the config_list is list(mut = list(), cnv = list(), exp = list(), ...). The list of each dataset must contain the following arguments
#'    list(dataset_name, 
#'         data_type_weight, 
#'         default_weight, 
#'         tumor_file, 
#'         cell_line_file)
#'  *  dataset_name: short name of comparison data type
#'  *  data_type_weight: a numeric weight for the data type (NOTE: data type weights must sum to 1)
#'  *  default_weight: default (background) weight for specified data type (EXP) (DEFAULT: 0.01)
#'  *  tumor_file: path to a tab delimited table file with HGNC gene symbol rownames and columns as tumour samples
#'  *  cell_line_file: path to a tab delimited table file with HGNC gene symbol rownames and columns as cell lines
#'  *  known_cancer_gene_weights_file: path to a file with weights for genes known
#'   to be recurrently altered in cancer (e.g. recurrently mutated genes in TCGA pan-cancer analyses). 
#'   A two-column tab-delimited file - the first column has the gene names and the second column specifies the weights (Default: NULL)
#'  *  cancer_specific_gene_weights_file: path to a file with weights for cancer-specific
#'   set of recurrently altered genes. A tab-delimited file - the first column has the gene names,
#'   and the second column specifies the weights (Default: NULL)
#' @param gene_list a vector of HGNC gene symbols to run comparison only for the specified genes (Default: NULL)
#' @param remove_errored_dataset_comparisons will skip the data types which 
#'   cannot be compared for technical reasons(no enaugh genes to compare, or 
#'   data contain only 0 values) when set to TRUE (Default: FALSE)
#' @param run_mds a boolean, whether to run multidimensional scaling (MDS) on dataset (Default: TRUE)
#' @param verbose show debugging information
#' 
#' @return a list with multiple items. Each 
#' * dist_mat: a matrix of combined pairwise distances for all data types
#' * dist_mat_by_data_type: a list of pairwise distances for each data type
#' * isomdsfit: a two-column (2-dimension) fitting of the distances reduced to two dimensions via MDS - multidimensional scaling, using the isoMDS function for all data types
#' * isomdsfit_by_data_type: a two-column (2-dimension) fitting of the distances reduced to two dimensions via MDS - multidimensional scaling, using the isoMDS function for each data type
#' * composite_mat: the composite matrix (see Details)
#' * cell_line_ids: a vector of cell line IDs/names with all data types
#' * tumor_ids: a vector of tumor IDs with all data types
#' * calculated_data_types: chacrachter vector of the dataset names which were analysed by comparison function
#' 
#' @examples 
#' tumor_mut_file <- system.file("extdata", "READ_data_for_running_TC", "tumor_mut.txt", 
#'   package="tumorcomparer")
#' tumor_cna_file <- system.file("extdata", "READ_data_for_running_TC", "tumor_cna.txt", 
#'   package="tumorcomparer")
#' tumor_exp_file <- system.file("extdata", "READ_data_for_running_TC", "tumor_exp.txt", 
#'   package="tumorcomparer")
#' 
#' cell_line_mut_file <- system.file("extdata", "READ_data_for_running_TC", "cell_line_mut.txt", 
#'   package="tumorcomparer")
#' cell_line_cna_file <- system.file("extdata", "READ_data_for_running_TC", "cell_line_cna.txt", 
#'   package="tumorcomparer")
#' cell_line_exp_file <- system.file("extdata", "READ_data_for_running_TC", "cell_line_exp.txt", 
#'   package="tumorcomparer")
#' 
#' known_cancer_gene_weights_mut_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
#' known_cancer_gene_weights_cna_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "default_weights_for_known_cancer_genes_cna.txt", package="tumorcomparer")
#' known_cancer_gene_weights_exp_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "default_weights_for_known_cancer_genes_exp.txt", package="tumorcomparer")
#' 
#' cancer_specific_gene_weights_mut_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "Genes_and_weights_mut.txt", package="tumorcomparer")
#' cancer_specific_gene_weights_cna_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "Genes_and_weights_cna.txt", package="tumorcomparer")
#' cancer_specific_gene_weights_exp_file <- system.file("extdata", "READ_data_for_running_TC", 
#'   "Genes_and_weights_exp.txt", package="tumorcomparer")
#'   
#' config_list <- list(
#'   mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
#'     tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
#'     known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
#'     cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
#'   cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
#'     tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
#'     known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
#'     cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
#'   exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
#'     tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
#'     known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
#'     cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
#' )
#' 
#' comparison_result <- run_comparison_config_list(config_list = config_list)
#' 
#' @export 
run_comparison_config_list <- function(config_list, 
                                       gene_list=NULL, 
                                       remove_errored_dataset_comparisons=FALSE, 
                                       run_mds=TRUE,
                                       verbose=FALSE) {
  
  ## Check input data
  if(length(config_list) < 1) {
    stop("ERROR: At least one data type must be provided")  
  }
  
  
  if(sum(sapply(config_list, function(x) {x$data_type_weight})) != 1) {
    stop("ERROR: Sum of *_data_type_weights must sum up to 1")  
  }
  
  ## Run comparison
  calculated_datasets_comparison <- lapply(config_list, function(x) {
    if(verbose) {
      cat("INFO: Data Type: ", x$dataset_name, "\n")      
    }

    tumor_dat <- read.table(x$tumor_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    cell_line_dat <- read.table(x$cell_line_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    
    if(sum(colSums(tumor_dat)) == 0 && sum(colSums(cell_line_dat)) == 0) {
      if(remove_errored_dataset_comparisons) {
        cat(paste0("INFO: Skipping ", x$dataset_name, " data; only zero values found", "\n"))
        
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
        gene_list=gene_list,
        verbose=verbose
        )
    }
  })
  
  ## filter skipped datasets
  calculated_datasets_comparison <- calculated_datasets_comparison[which(!sapply(calculated_datasets_comparison, is.null))]
  
  config_list <- config_list[names(calculated_datasets_comparison)]
  
  ## harmonize datasets
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
  if(run_mds) {
    isomdsfit <-  isoMDS(combined_dist, k=2)      
  } else {
    isomdsfit <- NA
  }
  
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
