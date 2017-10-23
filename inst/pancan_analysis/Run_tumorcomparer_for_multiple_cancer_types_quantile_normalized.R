require(tumorcomparer)

run_tc_qn_for_cancers_list <- function(cancers_list,gene_set_mut,gene_set_cna,gene_set_exp)
{
  # NECESSARY VARIABLES FROM Read_in_and_prepare.R ----
  # MUT: 
  mut_mat_TCGA_after_Annovar <- mut_mat_TCGA_after_Annovar
  mut_mat_CCLP_after_Annovar <- mut_mat_CCLP_after_Annovar
  
  ## CNA: 
  TCGA_Pancan_all_thresholded_by_genes_whitelisted <- TCGA_Pancan_all_thresholded_by_genes_whitelisted
  CCLP_GISTIC_thresholded <- CCLP_GISTIC_thresholded
  
  ## EXP for both CCLP and TCGA 
  CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized <- CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized
  
  ## Two columns from table (TBA)
  TCGA_id_and_tumor_type <- TCGA_id_and_tumor_type
  
  ##  RNASeq weights
  TCGA_PanCan_RNASeq_RPKMs_Weights <- TCGA_PanCan_RNASeq_RPKMs_Weights
  
  ## Assignment of cell line/tumor and cancer type 
  GDSC_1000_Binary_Event_Matrix <- GDSC_1000_Binary_Event_Matrix
  
  # SET BASIC PARAMETERS ----
  num_cancers <- length(cancers_list)
  #comparison_results_list <- vector('list', num_cancers)
  #categorization_results_list <- vector('list', num_cancers)
  comparison_results_list <- list()
  categorization_results_list <- list()
  
  # RUN ACROSS CANCER TYPES ----
  for(i in 1:num_cancers)
  {
    cancer_type <- cancers_list[i]
    
    # Create tumor mut, cna and exp files
    write.table(mut_mat_TCGA_after_Annovar[intersect(rownames(mut_mat_TCGA_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_TCGA_after_Annovar),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_mut.txt", sep = "\t", quote = FALSE)
    write.table(TCGA_Pancan_all_thresholded_by_genes_whitelisted[intersect(rownames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),gene_set_cna),intersect(colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_cna.txt", sep = "\t", quote = FALSE)
    write.table(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized[intersect(rownames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),gene_set_exp),intersect(colnames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_exp.txt", sep = "\t", quote = FALSE)
    
    # Create weight files
    write.table(file="Genes_and_weights_exp.txt",as.matrix(TCGA_PanCan_RNASeq_RPKMs_Weights)[,cancer_type],quote = FALSE,sep="\t")
    system(paste0("cp ",paste0("MUT_Weights_TCGA_",cancer_type,".txt") ," Genes_and_weights_mut.txt"))
    system(paste0("cp ",paste0("CNA_Weights_TCGA_",cancer_type,".txt") ," Genes_and_weights_cna.txt"))
    
    # Create cell line mut, cna and exp files
    if(cancer_type == "COAD" || cancer_type == "READ"){cancer_type <- "COAD/READ"} # CCLP uses COAD/READ
    if(cancer_type != "OV")
    {
      write.table(mut_mat_CCLP_after_Annovar[intersect(rownames(mut_mat_CCLP_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_CCLP_after_Annovar),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_mut.txt", sep = "\t", quote = FALSE)
      write.table(CCLP_GISTIC_thresholded[intersect(rownames(CCLP_GISTIC_thresholded),gene_set_cna),intersect(colnames(CCLP_GISTIC_thresholded),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_cna.txt", sep = "\t", quote = FALSE)
      write.table(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized[intersect(rownames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),gene_set_exp),intersect(colnames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_exp.txt", sep = "\t", quote = FALSE)
    }
    else
    {
      write.table(mut_mat_CCLP_after_Annovar[intersect(rownames(mut_mat_CCLP_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_CCLP_after_Annovar),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_mut.txt", sep = "\t", quote = FALSE)
      write.table(CCLP_GISTIC_thresholded[intersect(rownames(CCLP_GISTIC_thresholded),gene_set_cna),intersect(colnames(CCLP_GISTIC_thresholded),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_cna.txt", sep = "\t", quote = FALSE)
      write.table(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized[intersect(rownames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),gene_set_exp),intersect(colnames(CCLP_Microarray_TCGA_Tumors_RPKMs_quantile_normalized),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_exp.txt", sep = "\t", quote = FALSE)
    } 
    
    cat(cancer_type,":\n")
    system.time(source(system.file("pancan_analysis", "RunTumorComparerExample_updated_cont_exp.R", package="tumorcomparer")))
    comparison_results_list[[i]] <- comparison_result
    categorization_results_list[[i]] <- categorization_list
    names(comparison_results_list[[i]]) <- cancer_type
    names(categorization_results_list[[i]]) <- cancer_type
    
    #comparison_results_list <- list(comparison_results_list,comparison_result)
    #categorization_results_list <- list(categorization_results_list,categorization_list)
  }
  
  tc_results <- list(cancers_list = cancers_list, 
                     comparison_results_list_of_lists = comparison_results_list, 
                     categorization_results_list_of_lists = categorization_results_list)
  return(tc_results)
}
