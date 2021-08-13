# library(data.table)
# library(tumorcomparer)
# library(MASS)
# 
# ## checking compatible pathways for the genset comparisons
# gensets_and_available_genes <- as.data.frame(Reduce(rbind,
# 
#                                       lapply(names(tcga_pancan_pathway_genes), function(x) {
# 
#                                         t(
#                                           sapply(c("TCGA_Expression_Quantile_Normalized", "TCGA_GISTIC_all_data_by_genes", "mut_mat_TCGA_after_Annovar"), function(y) {
# 
#                                             cbind(x ,length(tcga_pancan_pathway_genes[[x]]), sum(tcga_pancan_pathway_genes[[x]] %in% get(y)$Gene))
# 
# 
#                                           }))
# 
#                                       })
# ), stringsAsFactors = F)
# 
# colnames(gensets_and_available_genes) <- c("geneset", "n_genes", "detected_genes")
# 
# gensets_and_available_genes$n_genes <- as.integer(gensets_and_available_genes$n_genes)
# gensets_and_available_genes$detected_genes <- as.integer(gensets_and_available_genes$detected_genes)
# 
# genesets_minimum_matches <- aggregate(gensets_and_available_genes$detected_genes, by = list(gensets_and_available_genes$geneset), min)
# 
# filtered_genesets <- unique(unname(genesets_minimum_matches$Group.1[which(genesets_minimum_matches$x > 4)]))
# 
# 
# ## precomputing all comparisons
# precomputed_comparisons <- lapply(as.character(avail_cancer_types), function(x) {
# 
#   lapply(filtered_genesets, function(y) {
# 
#     print(x)
#     print(y)
# 
#     comparison_list <- geneset_comparison(cancer_type = x, gene_list = tcga_pancan_pathway_genes[[y]], remove_errored_dataset_comparisons = TRUE)
# 
#     plot_data <- tumorcomparer::make_balloon_plot_data_from_comparison_result(comparison_list)
# 
#     list(compared_data_types = comparison_list$calculated_data_types, plot_data = plot_data)
# 
#   })
# 
# })
# 
# names(precomputed_comparisons) <- as.character(avail_cancer_types)
# 
# precomputed_comparisons <- lapply(precomputed_comparisons, function(x) {
# 
#   names(x) <- filtered_genesets
# 
#   x
# })
# 
# saveRDS(precomputed_comparisons, file = "../../data/precomputed_comparisons.rds")
# 
# ### detecting NaN cases
# 
# nan_datasets <- lapply(names(precomputed_comparisons), function(x) {
#   lapply(names(precomputed_comparisons[[x]]), function(y) {
#     if(is.nan(precomputed_comparisons[[x]][[y]]$plot_data[1,3])) {
#       c(x, y)
#     }
#   })
# })
# 
# nan_datasets <- Reduce(rbind, lapply(nan_datasets, function(x) {Reduce(rbind, x)}))
# 
# rownames(nan_datasets) <- NULL
# 
# colnames(nan_datasets) <- c("cancer", "geneset")
# 
# saveRDS(nan_datasets, file = "../temp_files/nan_datasets.rds")
# 
# 
# selected_geneset_comparions <- lapply(names(precomputed_comparisons), function(x) {
#   c("Most Variable Genes", 
#     unlist(lapply(names(precomputed_comparisons[[x]]), function(y) {
#       if(!is.nan(precomputed_comparisons[[x]][[y]]$plot_data[1,3])) {
#         y
#       }
#     }))
#     )
# })
# 
# names(selected_geneset_comparions) <- names(precomputed_comparisons)
# 
# saveRDS(selected_geneset_comparions, file = "data/selected_geneset_comparions.rds")
# 
# #### detected errors and special cases ####
# 
# ## mutation data with all zero values
# GBM_test <- geneset_comparison(cancer_type = "GBM", gene_list = tcga_pancan_pathway_genes[["MYC"]], remove_errored_dataset_comparisons = TRUE)
# 
# ## does not have matched cell lines
# geneset_comparison(cancer_type = "TGCT", gene_list = tcga_pancan_pathway_genes[["Cell Cycle"]], remove_errored_dataset_comparisons = TRUE)
# 
# 
# ## throws an error on cnv data while there are 5 genes present
# geneset_comparison(cancer_type = "LIHC", gene_list = tcga_pancan_pathway_genes[["TGF-Beta"]], remove_errored_dataset_comparisons = TRUE)
# 
