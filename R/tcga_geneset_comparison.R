#' geneset comparison function for TCGA dataset
#' 
#' @param tcga_dataset path to the tcga dataset (tcga_comparison_data.RData), see details
#' @param cancer_type TCGA cancer name (GBM, LGG, etc.)
#' @param gene_list a vector of HGNC gene symbols which will be used to perform the comparison
#' @param remove_errored_dataset_comparisons will skip the data types which can't be compared for technical reasons(no enaugh genes to compare, or data contain only 0 values) when set to TRUE (Default: FALSE)
#' 
#' @return returns comparison list for specified genes. The output list is described in run_comparison_config_list function documentation
#' 
#' @details
#' For running tcga_geneset_comparison you will need to download tcga_comparison_data.RData data file from ... and specify the path of the RData file in tcga_dataset argument
#' 
#' @importFrom data.table fwrite
#' 
#' @export 
tcga_geneset_comparison <- function(tcga_dataset, cancer_type, gene_list, remove_errored_dataset_comparisons = FALSE) {
  
  # FROM: https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  # to deal with "undefined global functions or variables" errors
  CCLP_Expression_Quantile_Normalized <- NULL
  CCLP_GISTIC_all_data_by_genes <- NULL
  CCLP_TCGA_Types_Combined <- NULL 
  Gene <- NULL
  TCGA_Expression_Quantile_Normalized <- NULL
  TCGA_GISTIC_all_data_by_genes <- NULL
  TCGA_id_and_tumor_type <- NULL
  mut_mat_CCLP_after_Annovar <- NULL
  mut_mat_TCGA_after_Annovar <- NULL
  
  load(tcga_dataset)
  
  if(length(cancer_type) > 1) {
    stop("ERROR: The analysis can be run only for one cancer type")  
  }
  
  most_variable_genes_precomputed_path <- system.file("extdata", 
                                                        "mtc_results_20200331", 
                                                        "mtc_results_20200331_no_factors.rds", 
                                                        package = "tumorcomparer")
  most_variable_genes_precomputed_results <- readRDS(most_variable_genes_precomputed_path)
  avail_cancer_types <- unique(most_variable_genes_precomputed_results$Tumor_Cancer_Type)
  
  if(!(cancer_type %in% avail_cancer_types)) {
    stop("ERROR: TCGA project with specified name is not found")
  }
  
  mut_cclp_ids <- c("Gene", intersect(colnames(mut_mat_CCLP_after_Annovar),CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)]))
  
  mut_tcga_ids <- c("Gene", intersect(colnames(mut_mat_TCGA_after_Annovar),TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)]))
  
  cna_cclp_ids <- c("Gene", intersect(colnames(CCLP_GISTIC_all_data_by_genes),CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)]))
  
  cna_tcga_ids <- c("Gene", intersect(colnames(TCGA_GISTIC_all_data_by_genes),TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)]))
  
  exp_cclp_ids <- c("Gene", intersect(colnames(CCLP_Expression_Quantile_Normalized),CCLP_TCGA_Types_Combined$Model_name[which(CCLP_TCGA_Types_Combined$TCGA_Type == cancer_type)]))
  
  exp_tcga_ids <- c("Gene", intersect(colnames(TCGA_Expression_Quantile_Normalized),TCGA_id_and_tumor_type$TCGA_ID[which(TCGA_id_and_tumor_type$Cancer_Type == cancer_type)]))
  
  if(sum(gene_list %in% Reduce(intersect, list(mut_mat_TCGA_after_Annovar$Gene, TCGA_GISTIC_all_data_by_genes$Gene, TCGA_Expression_Quantile_Normalized$Gene))) < 2) {
    
    stop("ERROR: The number of genes for comparison is less than 2")
    
  }
  
  ## data filtering
  fwrite(mut_mat_TCGA_after_Annovar[Gene %in% gene_list, ..mut_tcga_ids], file = "tumor_mut_filtered.txt", sep = "\t", quote = F)
  fwrite(TCGA_GISTIC_all_data_by_genes[Gene %in% gene_list, ..cna_tcga_ids], file = "tumor_cna_filtered.txt", sep = "\t", quote = F)
  fwrite(TCGA_Expression_Quantile_Normalized[Gene %in% gene_list, ..exp_tcga_ids], file = "tumor_exp_filtered.txt", sep = "\t", quote = F)
  
  
  fwrite(mut_mat_CCLP_after_Annovar[Gene %in% gene_list, ..mut_cclp_ids], file = "cell_line_mut_filtered.txt", sep = "\t", quote = F)
  fwrite(CCLP_GISTIC_all_data_by_genes[Gene %in% gene_list, ..cna_cclp_ids], file = "cell_line_cna_filtered.txt", sep = "\t", quote = F)
  fwrite(CCLP_Expression_Quantile_Normalized[Gene %in% gene_list, ..exp_cclp_ids], file = "cell_line_exp_filtered.txt", sep = "\t", quote = F)
  
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = "tumor_mut_filtered.txt", 
                               cell_line_file = "cell_line_mut_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_mut.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_mut.txt")
                               ),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = "tumor_cna_filtered.txt", 
                               cell_line_file = "cell_line_cna_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_cna.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_cna.txt")
                               ),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = "tumor_exp_filtered.txt", 
                               cell_line_file = "cell_line_exp_filtered.txt"
                               # known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_exp.txt"), 
                               # cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_exp.txt")
                               )
  )
  
  comparison_list <- run_comparison_config_list(config_list = config_list, remove_errored_dataset_comparisons = remove_errored_dataset_comparisons) 
  
  file.remove(c("tumor_mut_filtered.txt", "tumor_cna_filtered.txt", "tumor_exp_filtered.txt", 
                "cell_line_mut_filtered.txt", "cell_line_cna_filtered.txt", "cell_line_exp_filtered.txt"))
  
  return(comparison_list)
  
}
