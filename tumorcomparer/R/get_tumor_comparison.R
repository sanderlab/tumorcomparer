#' Run a comparison between 
get_tumor_comparison <- function(x) {
  # LOAD DATA ---- 
  CNA_default_weight <- 0.01
  MUT_default_weight <- 0.01
  
  CNA_known_cancer_gene_weight <- 0.1
  MUT_known_cancer_gene_weight <- 0.1
  
  tumor_mut_file <- "tumor_MUT.txt"
  tumor_cna_file <- "tumor_CNA.txt"
  cell_line_mut_file <- "cell_line_MUT.txt"
  cell_line_cna_file <- "cell_line_CNA.txt"
  
  pancancer_gene_weights_file <- "Default_weights_for_known_cancer_genes.txt"
  cancer_specific_gene_weights_file <- "Genes_and_weights.txt"
  
  output_composite_alteration_matrix_file <- "composite_alteration_matrix.txt"
  
  # GET INTERSECTING GENE ----
  tumor_MUT <- read.table(tumor_mut_file,sep="\t",header=TRUE,row.names=1)
  tumor_CNA <- read.table(tumor_cna_file,sep="\t",header=TRUE,row.names=1)
  cell_line_MUT <- read.table(cell_line_mut_file,sep="\t",header=TRUE,row.names=1)
  cell_line_CNA <- read.table(cell_line_cna_file,sep="\t",header=TRUE,row.names=1)
  
  tumors_with_both_MUT_and_CNA <- intersect(colnames(tumor_MUT),colnames(tumor_CNA))
  cell_lines_with_both_MUT_and_CNA <- intersect(colnames(cell_line_MUT),colnames(cell_line_CNA))
  cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, return_first_part)
  
  genes_with_MUT_in_both <- intersect(rownames(tumor_MUT),rownames(cell_line_MUT))
  genes_with_CNA_in_both <- intersect(rownames(tumor_CNA),rownames(cell_line_CNA))
  
  # Need not do this unless really want MUT and CNA data for the same gene sets
  genes_in_all_4_files <- intersect(genes_with_MUT_in_both,genes_with_CNA_in_both) 
  
  cell_line_CNA_high_level_only <- apply(cell_line_CNA,2,keep_only_high_level_cnas)
  tumor_CNA_high_level_only <- apply(tumor_CNA,2,keep_only_high_level_cnas)
  
  # CREATE COMPOSITE MATRIX ----
  ## Rows genes (MUT, CNA), columns (Tumors/cell lines)
  composite_CNA <- cbind(cell_line_CNA[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
  composite_CNA_high_level_only <- cbind(cell_line_CNA_high_level_only[genes_in_all_4_files,cell_lines_with_both_MUT_and_CNA],tumor_CNA_high_level_only[genes_in_all_4_files,tumors_with_both_MUT_and_CNA])
  composite_MUT <- cbind(cell_line_MUT[genes_with_MUT_in_both,cell_lines_with_both_MUT_and_CNA],tumor_MUT[genes_with_MUT_in_both,tumors_with_both_MUT_and_CNA])
  
  rownames(composite_MUT) <- paste(rownames(composite_MUT),"MUT",sep="_")
  rownames(composite_CNA) <- paste(rownames(composite_CNA),"CNA",sep="_")
  rownames(composite_CNA_high_level_only) <- paste(rownames(composite_CNA_high_level_only),"CNA",sep="_")
  
  composite_mat <- rbind(composite_MUT,composite_CNA)
  alt_mat <- rbind(composite_MUT,composite_CNA)
  
  # WRITE COMPOSITE
  write.table(composite_mat,file=output_composite_alteration_matrix_file,sep="\t",quote=FALSE)
  #composite_mat_high_level_only <- rbind(composite_MUT,composite_CNA_high_level_only)
  
  # Calculation of alteration frequencies 
  # Assign frequency weights as (freq. of alteration of gene)/(mean freq. of alteration across all genes) - "rewarding recurrent changes"
  overall_alt_freq <- length(which((alt_mat[]) != 0))/ ( length(which((alt_mat[])==0)) + length(which((alt_mat[])!=0)))
  
  freq_alt <- rep(0, nrow(alt_mat))
  freq_alt <- apply(alt_mat,1,compute_freq_alt)
  #freq_alt_high_level  <- apply(composite_mat_high_level_only,1,compute_freq_alt)
  freq_alt_mut_tumors <- apply(composite_MUT[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)
  freq_alt_cna_tumors <- apply(composite_CNA[,tumors_with_both_MUT_and_CNA],1,compute_freq_alt)
  
  freq_alt_samplewise <- apply(alt_mat,2,compute_freq_alt)
  #freq_alt_samplewise_CNA_high_level_only <- apply(composite_mat_high_level_only,2,compute_freq_alt)
  alt_mat <- alt_mat[,which(freq_alt_samplewise > 0)]
  composite_mat <- composite_mat[,which(freq_alt_samplewise > 0)]
  #composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise > 0)]
  
  names(freq_alt) <- rownames(alt_mat)
  freq_weights <- rep(1,nrow(alt_mat))
  freq_weights <- freq_alt/overall_alt_freq
  names(freq_weights) <- rownames(alt_mat)
  
  # GET WEIGHTS 
  # Read in user-provided weights
  known_cancer_genes_and_weights_all <- read.table(pancancer_gene_weights_file, sep="\t",header=TRUE,row.names=1) 
  known_cancer_genes_and_weights <- as.matrix(known_cancer_genes_and_weights_all[intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
  genes_and_weights_all <- read.table(cancer_specific_gene_weights_file,sep="\t",header=TRUE,row.names=1) 
  genes_and_weights <- as.matrix(genes_and_weights_all[intersect(rownames(genes_and_weights_all),rownames(alt_mat)),]) # To eliminate entries not present in alteration matrix, if any
  rownames(genes_and_weights) <- intersect(rownames(genes_and_weights_all),rownames(alt_mat))
  rownames(known_cancer_genes_and_weights) <- intersect(rownames(known_cancer_genes_and_weights_all),rownames(alt_mat))
  
  annotation_weights <- rep(NA,nrow(alt_mat)) 
  annotation_weights[grep("_CNA",rownames(alt_mat))] <- CNA_default_weight
  annotation_weights[grep("_MUT",rownames(alt_mat))] <- MUT_default_weight
  
  names(annotation_weights) <- rownames(alt_mat) 
  for(i in 1:nrow(known_cancer_genes_and_weights)) # Overwrite default weights for known cancer genes
    annotation_weights[rownames(known_cancer_genes_and_weights)[i]] = known_cancer_genes_and_weights[i,]
  for(i in 1:nrow(genes_and_weights)) # Overwrite default weights for input provided weights
    annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i,]
  
  
  gene_weights <- rep(1,nrow(alt_mat))
  names(gene_weights) <- rownames(alt_mat)
  
  gene_weights <- annotation_weights # if using user-provided weights only
  
  gene_weights <- gene_weights/max(gene_weights) # map to 0-1
  
  # CALCULATE CORRELATIONS ----
  # Including low-level CNAs
  cor_weighted <- calc_weighted_corr(as.matrix(composite_mat),as.matrix(composite_mat),gene_weights)
  # Excluding low-levels CNAs
  #cor_weighted_high_level_only <- calc_weighted_corr(as.matrix(composite_mat_high_level_only),as.matrix(composite_mat_high_level_only),gene_weights)
  cor_unweighted <- cor(alt_mat)

  results <- list(
    cor_weighted=cor_weighted, 
    cor_unweighted=cor_unweighted,
    composite_mat=composite_mat,
    cell_lines_with_both_MUT_and_CNA=cell_lines_with_both_MUT_and_CNA,
    tumors_with_both_MUT_and_CNA=tumors_with_both_MUT_and_CNA
  )
  
  return(results)
}