#' Gene Set Comparison Function for TCGA Dataset
#' 
#' @param gene_list a vector of HGNC gene symbols which will be used to perform the comparison
#' @param cancer_type TCGA cancer name (GBM, LGG, etc.)
#' @param tc_dataset path to RData for TumorComparer dataset from publication; see details
#' @param tcga_ids a vector of TCGA IDs for analysis; cannot be combined with cancer_type
#' @param cclp_ids a vector of CCLP IDs for analysis; cannot be combined with cancer_type
#' @param tc_dataset_dir path to directory with TumorComparer dataset as text files (recommended instead of tc_dataset); see details
#' @param remove_errored_dataset_comparisons will skip the data types which can't be compared for technical reasons(no enaugh genes to compare, or data contain only 0 values) when set to TRUE (Default: FALSE)
#' @param remove_tmp_files remove temporary files (Default: TRUE)
#' @param run_mds a boolean, whether to run multidimensional scaling (MDS) on dataset (Default: TRUE)
#' @param verbose show debugging information
#' 
#' @return returns comparison list for specified genes. The output list is described in run_comparison_config_list function documentation
#' 
#' @note # FIXME ADD EXAMPLE
#' 
#' @details
#' This function requires downloading: https://zenodo.org/record/4627644/files/tc_data_pancancer_march2021.tar.gz?download=1
#' 
#' The above data is the dataset used in the 2021 TumorComparer publication.  
#' Data can then be loaded from the extracted text files if 'tc_dataset_dir' 
#' parameter is provided or data can be loaded from an RData workspace variable
#' using the 'tc_dataset' parameter if it has been pre-generated. tc_dataset_dir' 
#' is recommended.
#' 
#' @importFrom data.table fwrite fread
#' 
#' @export 
tc_geneset_comparison <- function(gene_list, 
                                  cancer_type=NULL, 
                                  tcga_ids=NULL,
                                  cclp_ids=NULL, 
                                  tc_dataset=NULL, 
                                  tc_dataset_dir=NULL,
                                  remove_errored_dataset_comparisons=FALSE, 
                                  remove_tmp_files=TRUE,
                                  run_mds=TRUE,
                                  verbose=FALSE) {
  
  # FROM: https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  # to deal with "undefined global functions or variables" errors
  CCLP_Expression_Quantile_Normalized <- NULL
  CCLP_GISTIC_all_data_by_genes <- NULL
  CCLP_TCGA_Types_Combined <- NULL 
  TCGA_Expression_Quantile_Normalized <- NULL
  TCGA_GISTIC_all_data_by_genes <- NULL
  TCGA_id_and_tumor_type <- NULL
  mut_mat_CCLP_after_Annovar <- NULL
  mut_mat_TCGA_after_Annovar <- NULL
  Gene <- NULL
  
  if(!is.null(tc_dataset_dir)) {
    TCGA_Expression_Quantile_Normalized <- fread(input=file.path(tc_dataset_dir, "TCGA_Expression_Quantile_Normalized.txt"))
    colnames(TCGA_Expression_Quantile_Normalized)[1] <- "Gene"
    
    TCGA_GISTIC_all_data_by_genes <- fread(input=file.path(tc_dataset_dir, "TCGA_GISTIC_all_data_by_genes.txt"))
    colnames(TCGA_GISTIC_all_data_by_genes)[1] <- "Gene"
    
    CCLP_Expression_Quantile_Normalized <- fread(input=file.path(tc_dataset_dir, "CCLP_Expression_Quantile_Normalized.txt"))
    colnames(CCLP_Expression_Quantile_Normalized)[1] <- "Gene"
    
    CCLP_GISTIC_all_data_by_genes <- fread(input=file.path(tc_dataset_dir, "CCLP_GISTIC_all_data_by_genes.txt"))
    colnames(CCLP_GISTIC_all_data_by_genes)[1] <- "Gene"
    
    mut_mat_CCLP_after_Annovar <- fread(input=file.path(tc_dataset_dir, "CCLP_mutation_matrix.txt"))
    colnames(mut_mat_CCLP_after_Annovar)[1] <- "Gene"
    
    mut_mat_TCGA_after_Annovar <- fread(input=file.path(tc_dataset_dir, "TCGA_mutation_matrix.txt"))
    colnames(mut_mat_TCGA_after_Annovar)[1] <- "Gene"
    
    CCLP_TCGA_Types_Combined <- fread(input=file.path(tc_dataset_dir, "CCLP_ID_and_Cancer_Type.txt"))
    TCGA_id_and_tumor_type <- fread(input=file.path(tc_dataset_dir, "TCGA_ID_and_Cancer_Type.txt"))
  } else if(!is.null(tc_dataset)) {
    load(tc_dataset)
  } else {
    stop("ERROR: TumorComparer dataset location was not specified")
  }
  
  # TODO: Create test for this conditions
  if(is.null(cancer_type) && !is.null(tcga_ids) && !is.null(cclp_ids)) {
    tmp_mut_cclp_ids <- cclp_ids
    tmp_mut_tcga_ids <- tcga_ids
    tmp_cna_cclp_ids <- cclp_ids
    tmp_cna_tcga_ids <- tcga_ids
    tmp_exp_cclp_ids <- cclp_ids
    tmp_exp_tcga_ids <- tcga_ids
  } else if(!is.null(cancer_type) && is.null(tcga_ids) && is.null(cclp_ids)) {
    if(length(cancer_type) > 1) {
      stop("ERROR: The analysis can be run only for one cancer type")  
    }
    
    most_variable_genes_precomputed_path <- system.file("extdata", 
                                                        "mtc_results_20200331", 
                                                        "mtc_results_20200331_no_factors.rds", 
                                                        package = "tumorcomparer")
    most_variable_genes_precomputed_results <- readRDS(most_variable_genes_precomputed_path)
    avail_cancer_types <- as.character(unique(most_variable_genes_precomputed_results$Tumor_Cancer_Type))
    
    if(!(cancer_type %in% avail_cancer_types)) {
      stop("ERROR: TCGA cancer type with specified name is not found")
    }
    
    tmp_mut_cclp_ids <- intersect(colnames(mut_mat_CCLP_after_Annovar), CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)])
    tmp_mut_tcga_ids <- intersect(colnames(mut_mat_TCGA_after_Annovar), TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)])
    tmp_cna_cclp_ids <- intersect(colnames(CCLP_GISTIC_all_data_by_genes), CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)])
    tmp_cna_tcga_ids <- intersect(colnames(TCGA_GISTIC_all_data_by_genes), TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)])
    tmp_exp_cclp_ids <- intersect(colnames(CCLP_Expression_Quantile_Normalized), CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)])
    tmp_exp_tcga_ids <- intersect(colnames(TCGA_Expression_Quantile_Normalized), TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)])
  } else {
    stop("ERROR: Either (cancer_type) OR (tcga_ids and cclp_ids) must be set, but not both.")
  }
  
  mut_cclp_ids <- c("Gene", tmp_mut_cclp_ids)
  mut_tcga_ids <- c("Gene", tmp_mut_tcga_ids)
  cna_cclp_ids <- c("Gene", tmp_cna_cclp_ids)
  cna_tcga_ids <- c("Gene", tmp_cna_tcga_ids)
  exp_cclp_ids <- c("Gene", tmp_exp_cclp_ids)
  exp_tcga_ids <- c("Gene", tmp_exp_tcga_ids)
  
  if(sum(gene_list %in% Reduce(intersect, list(mut_mat_TCGA_after_Annovar$Gene, TCGA_GISTIC_all_data_by_genes$Gene, TCGA_Expression_Quantile_Normalized$Gene))) < 2) {
    stop("ERROR: The number of genes for comparison is less than 2")
  }
  
  ## Filter data
  fwrite(mut_mat_TCGA_after_Annovar[Gene %in% gene_list, mut_tcga_ids, with = FALSE], file = "tumor_mut_filtered.txt", sep = "\t", quote = FALSE)
  fwrite(TCGA_GISTIC_all_data_by_genes[Gene %in% gene_list, cna_tcga_ids, with = FALSE], file = "tumor_cna_filtered.txt", sep = "\t", quote = FALSE)
  fwrite(TCGA_Expression_Quantile_Normalized[Gene %in% gene_list, exp_tcga_ids, with = FALSE], file = "tumor_exp_filtered.txt", sep = "\t", quote = FALSE)
  
  fwrite(mut_mat_CCLP_after_Annovar[Gene %in% gene_list, mut_cclp_ids, with = FALSE], file = "cell_line_mut_filtered.txt", sep = "\t", quote = FALSE)
  fwrite(CCLP_GISTIC_all_data_by_genes[Gene %in% gene_list, cna_cclp_ids, with = FALSE], file = "cell_line_cna_filtered.txt", sep = "\t", quote = FALSE)
  fwrite(CCLP_Expression_Quantile_Normalized[Gene %in% gene_list, exp_cclp_ids, with = FALSE], file = "cell_line_exp_filtered.txt", sep = "\t", quote = FALSE)
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 1, 
                               tumor_file = "tumor_mut_filtered.txt", 
                               cell_line_file = "cell_line_mut_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_mut.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_mut.txt")
                               ),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 1, 
                               tumor_file = "tumor_cna_filtered.txt", 
                               cell_line_file = "cell_line_cna_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_cna.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_cna.txt")
                               ),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 1, 
                               tumor_file = "tumor_exp_filtered.txt", 
                               cell_line_file = "cell_line_exp_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_exp.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_exp.txt")
                               )
  )
  
  comparison_list <- run_comparison_config_list(config_list=config_list, 
                                                remove_errored_dataset_comparisons=remove_errored_dataset_comparisons, 
                                                run_mds=run_mds,
                                                verbose=verbose) 
  
  if(remove_tmp_files) {
    rm(TCGA_Expression_Quantile_Normalized,
       TCGA_GISTIC_all_data_by_genes,
       CCLP_Expression_Quantile_Normalized, 
       CCLP_GISTIC_all_data_by_genes, 
       mut_mat_CCLP_after_Annovar, 
       mut_mat_TCGA_after_Annovar)
    
    file.remove(c("tumor_mut_filtered.txt", "tumor_cna_filtered.txt", 
                  "tumor_exp_filtered.txt", "cell_line_mut_filtered.txt", 
                  "cell_line_cna_filtered.txt", "cell_line_exp_filtered.txt"))    
  }

  return(comparison_list)
}
