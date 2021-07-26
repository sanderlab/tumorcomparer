library(data.table)
library(tumorcomparer)
library(MASS)

## checking compatible pathways for the genset comparisons
gensets_and_available_genes <- as.data.frame(Reduce(rbind, 
                                      
                                      lapply(names(tcga_pancan_pathway_genes), function(x) {
                                        
                                        t(
                                          sapply(c("TCGA_Expression_Quantile_Normalized", "TCGA_GISTIC_all_data_by_genes", "mut_mat_TCGA_after_Annovar"), function(y) {
                                            
                                            cbind(x ,length(tcga_pancan_pathway_genes[[x]]), sum(tcga_pancan_pathway_genes[[x]] %in% get(y)$Gene))
                                            
                                            
                                          }))
                                        
                                      })
), stringsAsFactors = F)

colnames(gensets_and_available_genes) <- c("geneset", "n_genes", "detected_genes")

gensets_and_available_genes$n_genes <- as.integer(gensets_and_available_genes$n_genes)
gensets_and_available_genes$detected_genes <- as.integer(gensets_and_available_genes$detected_genes)

genesets_minimum_matches <- aggregate(gensets_and_available_genes$detected_genes, by = list(gensets_and_available_genes$geneset), min)

filtered_genesets <- unique(unname(genesets_minimum_matches$Group.1[which(genesets_minimum_matches$x > 4)]))


## precomputing all comparisons
precomputed_comparisons <- lapply(unique(CCLP_TCGA_Types_Combined$TCGA_Type), function(x) {
  
  lapply(filtered_genesets, function(y) {
    
    print(x)
    print(y)
    
    comparison_list <- geneset_comparison(cancer_type = x, gene_list = tcga_pancan_pathway_genes[[y]], remove_errored_dataset_comparisons = TRUE)
    
    plot_data <- tumorcomparer::make_balloon_plot_data_from_comparison_result(comparison_list)
    
    list(compared_data_types = comparison_list$calculated_data_types, plot_data = plot_data)
    
  })
  
})

names(precomputed_comparisons) <- unique(TCGA_id_and_tumor_type$Cancer_Type)

saveRDS(precomputed_comparisons, file = "../../data/precomputed_comparisons.rds")



#### detected errors and special cases ####

## mutation data with all zero values
geneset_comparison(cancer_type = "GBM", gene_list = tcga_pancan_pathway_genes[["MYC"]])

## does not have matched cell lines
geneset_comparison(cancer_type = "TGCT", gene_list = tcga_pancan_pathway_genes[["Cell Cycle"]], remove_errored_dataset_comparisons = TRUE)


## throws an error on cnv data while there are 5 genes present
geneset_comparison(cancer_type = "LIHC", gene_list = tcga_pancan_pathway_genes[["TGF-Beta"]], remove_errored_dataset_comparisons = TRUE)

