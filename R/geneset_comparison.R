library(data.table)

# TCGA ID and tumor types from TCGA PanCancer annotations
TCGA_id_and_tumor_type <- fread("../TC_Data_PanCancer_March2021/TCGA_ID_and_Cancer_Type.txt",check.names=F, nThread = 4)
# CCLP ID and TCGA tumor type from Cell Model Passport annotations
CCLP_TCGA_Types_Combined  <- fread("../TC_Data_PanCancer_March2021/CCLP_ID_and_Cancer_Type.txt",check.names=F, nThread = 4)

## reading mutation data
mut_mat_TCGA_after_Annovar <- fread("../TC_Data_PanCancer_March2021/TCGA_mutation_matrix.txt",check.names=F, nThread = 4)
mut_mat_CCLP_after_Annovar <- fread("../TC_Data_PanCancer_March2021/CCLP_mutation_matrix.txt",check.names=F, nThread = 4)

## unifiing cell line and TCGA genes
mut_mat_TCGA_after_Annovar <- mut_mat_TCGA_after_Annovar[Gene %in% intersect(mut_mat_TCGA_after_Annovar$Gene, mut_mat_CCLP_after_Annovar$Gene)]
mut_mat_CCLP_after_Annovar <- mut_mat_CCLP_after_Annovar[Gene %in% intersect(mut_mat_TCGA_after_Annovar$Gene, mut_mat_CCLP_after_Annovar$Gene)]

## reading CNV data
TCGA_GISTIC_all_data_by_genes <- fread("../TC_Data_PanCancer_March2021/TCGA_GISTIC_all_data_by_genes.txt",check.names=F, nThread = 4)
colnames(TCGA_GISTIC_all_data_by_genes)[1] <- "Gene"
CCLP_GISTIC_all_data_by_genes <- fread("../TC_Data_PanCancer_March2021/CCLP_GISTIC_all_data_by_genes.txt",check.names=F, nThread = 4)
colnames(CCLP_GISTIC_all_data_by_genes)[1] <- "Gene"

## unifiing cell line and TCGA genes
TCGA_GISTIC_all_data_by_genes <- TCGA_GISTIC_all_data_by_genes[Gene %in% intersect(TCGA_GISTIC_all_data_by_genes$Gene, CCLP_GISTIC_all_data_by_genes$Gene)]
CCLP_GISTIC_all_data_by_genes <- CCLP_GISTIC_all_data_by_genes[Gene %in% intersect(TCGA_GISTIC_all_data_by_genes$Gene, CCLP_GISTIC_all_data_by_genes$Gene)]

## reading expression data
TCGA_Expression_Quantile_Normalized <- fread("../TC_Data_PanCancer_March2021/TCGA_Expression_Quantile_Normalized.txt",check.names=F, nThread = 4)
colnames(TCGA_Expression_Quantile_Normalized)[1] <- "Gene"
CCLP_Expression_Quantile_Normalized <- fread("../TC_Data_PanCancer_March2021/CCLP_Expression_Quantile_Normalized.txt",check.names=F, nThread = 4)

## unification is not needed as the dat is alredy unified

geneset_comparison <- function(cancer_type, gene_list) {
  
  if(length(cancer_type) > 1) {
    
    stop("ERROR: The analysis can be run only for one cancer type")  
    
  }
  
  if(!(cancer_type %in% unique(TCGA_id_and_tumor_type$Cancer_Type))) {
    
    stop("ERROR: TCGA projec with specified name is not fount")  
    
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
                               cell_line_file = "cell_line_mut_filtered.txt",
                               known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_mut.txt"), 
                               cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_mut.txt")),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = "tumor_cna_filtered.txt", 
                               cell_line_file = "cell_line_cna_filtered.txt",
                               known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_cna.txt"), 
                               cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_cna.txt")),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = "tumor_exp_filtered.txt", 
                               cell_line_file = "cell_line_exp_filtered.txt",
                               known_cancer_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/default_weights_for_known_cancer_genes_exp.txt"), 
                               cancer_specific_gene_weights_file = paste0("../TC_Data_For_Each_Cancer_Type_March2021/", cancer_type, "/genes_and_weights_exp.txt"))
  )
  
  comparison_list <- run_comparison_config_list(config_list = config_list) 
  
  file.remove(c("tumor_mut_filtered.txt", "tumor_cna_filtered.txt", "tumor_exp_filtered.txt", 
                "cell_line_mut_filtered.txt", "cell_line_cna_filtered.txt", "cell_line_exp_filtered.txt"))
  
  return(comparison_list)
  
}


test_comparison <- geneset_comparison(cancer_type = "LGG", gene_list = tcga_pancan_pathway_genes$HIPPO)

