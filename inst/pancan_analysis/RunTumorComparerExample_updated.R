# Run example 
#setwd("~/Google Drive/example_files")

comparison_result <- run_comparison(available_data_types=c("mut", "cna", "exp"), 
               mut_data_type_weight = 1/3,
               cna_data_type_weight = 1/3,
               exp_data_type_weight = 1/3,
               cna_default_weight=0.01, 
               mut_default_weight=0.01,
               exp_default_weight=0.01,
               cna_known_cancer_gene_weight=0.1, 
               mut_known_cancer_gene_weight=0.1, 
               exp_known_cancer_gene_weight=0.1, 
               tumor_mut_file="tumor_mut.txt", 
               tumor_cna_file="tumor_cna.txt", 
               tumor_exp_file="tumor_exp.txt", 
               cell_line_mut_file="cell_line_mut.txt", 
               cell_line_cna_file="cell_line_cna.txt", 
               cell_line_exp_file="cell_line_exp.txt", 
               known_cancer_gene_weights_mut_file="default_weights_for_known_cancer_genes_mut.txt", 
               known_cancer_gene_weights_cna_file="default_weights_for_known_cancer_genes_cna.txt", 
               known_cancer_gene_weights_exp_file="default_weights_for_known_cancer_genes_exp.txt", 
               cancer_specific_gene_weights_mut_file="Genes_and_weights_mut.txt", 
               cancer_specific_gene_weights_cna_file="Genes_and_weights_cna.txt", 
               cancer_specific_gene_weights_exp_file="Genes_and_weights_exp.txt",
               #distance_similarity_measures=c("generalized_jaccard", "generalized_jaccard", "weighted_correlation"))
               distance_similarity_measures=c("generalized_jaccard", "generalized_jaccard", "generalized_jaccard"))

categorization_list <- categorize_cell_lines(num_tumors_for_comparison= length(comparison_result$tumor_ids) - 1,
                                  comparison_result$dist_mat,
                                  comparison_result$cell_line_ids,
                                  comparison_result$tumor_ids,
                                  trim_cell_line_names=FALSE) 

plot_mds(comparison_result,
                     categorization_list,
                     trim_cell_line_names=FALSE,
                     tumor_color="blue",
                     cell_line_color="orange",
                     use_gradient=TRUE,
                     tumor_shape=20,
                     cell_line_shape=17)
