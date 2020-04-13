run_tc_for_cancers_list <- function(cancers_list,gene_set_mut,gene_set_cna,gene_set_exp)
{
  num_cancers <- length(cancers_list)
  #comparison_results_list <- vector('list', num_cancers)
  #categorization_results_list <- vector('list', num_cancers)
  comparison_results_list <- list()
  categorization_results_list <- list()
  
  for(i in 1:num_cancers)
  {
    cancer_type <- cancers_list[i]
    
    # Create tumor mut, cna and exp files
    write.table(mut_mat_TCGA_after_Annovar[intersect(rownames(mut_mat_TCGA_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_TCGA_after_Annovar),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_mut.txt", sep = "\t", quote = F)
    write.table(TCGA_Pancan_all_thresholded_by_genes_whitelisted[intersect(rownames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),gene_set_cna),intersect(colnames(TCGA_Pancan_all_thresholded_by_genes_whitelisted),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_cna.txt", sep = "\t", quote = F)
    write.table(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores[intersect(rownames(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores),gene_set_exp),intersect(colnames(TCGA_PanCan_RNASeq_RPKMs_logged_Tumors_Discretized_By_Zscores),TCGA_id_and_tumor_type[which(TCGA_id_and_tumor_type[,2] == cancer_type),1])],file="tumor_exp.txt", sep = "\t", quote = F)
    
    # Create weight files
    write.table(file="Genes_and_weights_exp.txt",as.matrix(TCGA_PanCan_RNASeq_RPKMs_Weights)[,cancer_type],quote = F,sep="\t")
    system(paste0("cp ",paste0("MUT_Weights_TCGA_",cancer_type,".txt") ," Genes_and_weights_mut.txt"))
    system(paste0("cp ",paste0("CNA_Weights_TCGA_",cancer_type,".txt") ," Genes_and_weights_cna.txt"))
    
    # Create cell line mut, cna and exp files
    if(cancer_type == "COAD" || cancer_type == "READ"){cancer_type <- "COAD/READ"} # CCLP uses COAD/READ
      if(cancer_type != "OV")
      {
        write.table(mut_mat_CCLP_after_Annovar[intersect(rownames(mut_mat_CCLP_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_CCLP_after_Annovar),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_mut.txt", sep = "\t", quote = F)
        write.table(CCLP_GISTIC_thresholded[intersect(rownames(CCLP_GISTIC_thresholded),gene_set_cna),intersect(colnames(CCLP_GISTIC_thresholded),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_cna.txt", sep = "\t", quote = F)
        write.table(CCLP_Expression_Microarray_Discretized_By_Zscores[intersect(rownames(CCLP_Expression_Microarray_Discretized_By_Zscores),gene_set_exp),intersect(colnames(CCLP_Expression_Microarray_Discretized_By_Zscores),colnames(GDSC_1000_Binary_Event_Matrix)[which(GDSC_1000_Binary_Event_Matrix[1,] == cancer_type & GDSC_1000_Binary_Event_Matrix[2,] == "cell line")])],file="cell_line_exp.txt", sep = "\t", quote = F)
      }
      else
      {
        write.table(mut_mat_CCLP_after_Annovar[intersect(rownames(mut_mat_CCLP_after_Annovar),gene_set_mut),intersect(colnames(mut_mat_CCLP_after_Annovar),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_mut.txt", sep = "\t", quote = F)
        write.table(CCLP_GISTIC_thresholded[intersect(rownames(CCLP_GISTIC_thresholded),gene_set_cna),intersect(colnames(CCLP_GISTIC_thresholded),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_cna.txt", sep = "\t", quote = F)
        write.table(CCLP_Expression_Microarray_Discretized_By_Zscores[intersect(rownames(CCLP_Expression_Microarray_Discretized_By_Zscores),gene_set_exp),intersect(colnames(CCLP_Expression_Microarray_Discretized_By_Zscores),CCLP_Cell_Lines_And_Tissue_Types$`Cell lines`[which(CCLP_Cell_Lines_And_Tissue_Types$`Tissue type` == "ovary")])],file="cell_line_exp.txt", sep = "\t", quote = F)
      } 
      
      cat(cancer_type,":\n")
      system.time(source(system.file("pancan_analysis", "RunTumorComparerExample_updated.R", package="tumorcomparer")))
      comparison_results_list[[i]] <- comparison_result
      categorization_results_list[[i]] <- categorization_list
      names(comparison_results_list[[i]]) <- cancer_type
      names(categorization_results_list[[i]]) <- cancer_type
      
      #comparison_results_list <- list(comparison_results_list,comparison_result)
      #categorization_results_list <- list(categorization_results_list,categorization_list)
  }
  
  tc_results <- list(cancers_list = cancers_list,comparison_results_list_of_lists = comparison_results_list, categorization_results_list_of_lists = categorization_results_list)
  
  return(tc_results)
}
