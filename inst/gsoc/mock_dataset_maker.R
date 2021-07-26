set.seed(1)

selected_genes <- sample(sample(rownames(read.table(system.file("extdata", "ovarian_tcga_cclp", "tumor_mut.txt", package="tumorcomparer"), sep = "\t")), size = 100), size = 100)
  
selected_tumor_samples <- sample(colnames(read.table(system.file("extdata", "ovarian_tcga_cclp", "tumor_mut.txt", package="tumorcomparer"), sep = "\t")), size = 100)
  
selected_cell_line_samples <- colnames(read.delim(system.file("extdata", "ovarian_tcga_cclp", "cell_line_mut.txt", package="tumorcomparer"), sep = "\t"))
  
dataset_list <- c("exp", "cna", "meth", "mut", "prot")
  
mock_dataset <- lapply(dataset_list, function(x) {
    
  if(x %in% c("mut", "meth")) {
      
    tumor_file <- matrix(sample(c(1,0), size = 10000, replace = TRUE), 100, 100)
      
    cell_line_file <- matrix(sample(c(1,0), size = 100*length(selected_cell_line_samples), replace = TRUE), 100, length(selected_cell_line_samples))
      
  } else {
    
    tumor_file <- matrix(abs(rnorm(10000) + 10)*10, 100, 100)
    
    cell_line_file <- matrix(sample(c(1,0), size = 100*length(selected_cell_line_samples), replace = TRUE), 100, length(selected_cell_line_samples))
  
  }
    
  rownames(tumor_file) <- selected_genes
  rownames(cell_line_file) <- selected_genes
  
  colnames(tumor_file) <- selected_tumor_samples
  colnames(cell_line_file) <- selected_cell_line_samples
    
    
  write.table(tumor_file, file = paste0("../extdata/mock_5_data_types/tumor_", x, ".txt"), sep = "\t", row.names = T, col.names = NA)
    
  write.table(cell_line_file, file = paste0("../extdata/mock_5_data_types/cell_line_", x, ".txt"), sep = "\t", row.names = T, col.names = NA)
    
  
  known_cancer_gene_weights_file <- read.table(system.file("extdata", "ovarian_tcga_cclp", "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    
  cancer_specific_gene_weights_file <- read.table(system.file("extdata", "ovarian_tcga_cclp", "Genes_and_weights_mut.txt", package="tumorcomparer"), sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    
  write.table(known_cancer_gene_weights_file, file = paste0("../extdata/mock_5_data_types/default_weights_for_known_cancer_genes_", x, ".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
    
  write.table(cancer_specific_gene_weights_file, file = paste0("../extdata/mock_5_data_types/Genes_and_weights_", x, ".txt"), sep = "\t", row.names = T, col.names = NA, quote = F)
  
})