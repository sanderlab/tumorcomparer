#' Run a comparison between
#'
#' @param CNA_default_weight default weight for copy number alterations (CNA) (DEFAULT: 0.01)
#' @param MUT_default_weight default weight for mutation alterations (MUT) (DEFAULT: 0.01)
#' @param CNA_known_cancer_gene_weight a default weight (DEFAULT: 0.1) - genes in the TCGA pan-cancer CNA list (GISTIC peaks)
#' @param MUT_known_cancer_gene_weight a default weight (DEFAULT: 0.1) - genes in the TCGA pan-cancer mutation list (MUTSIG SMGs)
#' @param tumor_mut_file a file with binary mutation data for tumors, over the 1651 genes profiled by CCLE
#' @param tumor_cna_file a file with 5-valued GISTIC data for tumors, over 1529 genes (subset of 1651 genes above)
#' @param cell_line_mut_file a file with binary mutation data for cell lines, over the 1651 genes profiled by CCLE
#' @param cell_line_cna_file a file with 5-valued GISTIC data for cell lines, over 1529 genes (subset of 1651 genes above)
#' @param pancancer_gene_weights_file - a file with weights for the TCGA pan-cancer set of recurrent mutations and CNAs (copy number alterations)
#' @param cancer_specific_gene_weights_file a file with weights for cancer-specific set of recurrent mutations and CNAs (copy number alterations)
#' @param output_composite_alteration_matrix_file a string with filename for output of the the composite alteration matrix (see Details) (DEFAULT: "composite_alteration_matrix.txt") 
#' @param distance_similarity_measure (See Details) 
#'   (OPTIONS: "weighted_correlation", "generalized_jaccard") 
#' 
#' @details The composite matrix is a single matrix where the columns are samples 
#'   (i.e. tumors AND cell line IDs) and the rows are an rbind() of mutations 
#'   (with 1 or 0 outputs for each sample) and copy number alterations from 
#'   GISTIC (with values -2, -1, 0, 1, 2). Available similarity/distance measures include: 
#'   \itemize{
#'   \item{"weighted_correlation"}{Weighted correlation, based on weighted means and standard deviations}
#'   \item{"generalized_jaccard"}{A weighted distance based on the Jaccard coefficient}
#'    }
#'    
#' @return a list with multiple items. NOTE: The values of the dist and isomdsfit will depend on parameter "distance_similarity_measure".
#' \itemize{
#'   \item{"dist"}{a matrix of pairwise distances}
#'   \item{"isomdsfit"}{a two-column (2-dimension) fitting of the distances reduced to two dimensions via MDS - multidimensional scaling}
#'   \item{"cor_unweighted"}{a matrix of unweighted pairwise correlations}
#'   \item{"composite_mat"}{the composite matrix (see Details)}
#'   \item{"cell_lines_with_both_MUT_and_CNA"}{a vector of cell lines with both mutation (MUT) and copy number alteration (CNA) information}
#'   \item{"tumors_with_both_MUT_and_CNA"}{a vector of cell lines with both mutation (MUT) and copy number alteration (CNA) information}
#' }
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
get_tumor_comparison <- function(CNA_default_weight=0.01, 
                                 MUT_default_weight=0.01,
                                 NA_known_cancer_gene_weight=0.1, 
                                 MUT_known_cancer_gene_weight=0.1, 
                                 tumor_mut_file="tumor_MUT.txt", 
                                 tumor_cna_file="tumor_CNA.txt", 
                                 cell_line_mut_file="cell_line_MUT.txt", 
                                 cell_line_cna_file="cell_line_CNA.txt", 
                                 pancancer_gene_weights_file="Default_weights_for_known_cancer_genes.txt", 
                                 cancer_specific_gene_weights_file="Genes_and_weights.txt", 
                                 output_composite_alteration_matrix_file="composite_alteration_matrix.txt",
                                 distance_similarity_measure=c("weighted_correlation", "generalized_jaccard")) {

  # GET INTERSECTING GENES BETWEEN TUMORS AND CELL LINES ----
  tumor_MUT <- read.table(tumor_mut_file, sep = "\t", header = TRUE, row.names = 1)
  tumor_CNA <- read.table(tumor_cna_file, sep = "\t", header = TRUE, row.names = 1)
  cell_line_MUT <- read.table(cell_line_mut_file, sep = "\t", header = TRUE, row.names = 1)
  cell_line_CNA <- read.table(cell_line_cna_file, sep = "\t", header = TRUE, row.names = 1)
  
  tumors_with_both_MUT_and_CNA <- intersect(colnames(tumor_MUT), colnames(tumor_CNA))
  cell_lines_with_both_MUT_and_CNA <- intersect(colnames(cell_line_MUT), colnames(cell_line_CNA))
  cell_line_ids <- sapply(cell_lines_with_both_MUT_and_CNA, return_first_part)
  
  genes_with_MUT_in_both <- intersect(rownames(tumor_MUT), rownames(cell_line_MUT))
  genes_with_CNA_in_both <- intersect(rownames(tumor_CNA), rownames(cell_line_CNA))
  
  # Need not do this unless really want MUT and CNA data for the same gene sets
  genes_in_all_4_files <-
    intersect(genes_with_MUT_in_both, genes_with_CNA_in_both)
  
  cell_line_CNA_high_level_only <- apply(cell_line_CNA, 2, keep_only_high_level_cnas)
  tumor_CNA_high_level_only <- apply(tumor_CNA, 2, keep_only_high_level_cnas)
  
  # CREATE COMPOSITE MATRIX ----
  ## Rows genes (MUT, CNA), columns (Tumors/cell lines)
  composite_CNA <- cbind(cell_line_CNA[genes_in_all_4_files, cell_lines_with_both_MUT_and_CNA], tumor_CNA[genes_in_all_4_files, tumors_with_both_MUT_and_CNA])
  composite_CNA_high_level_only <- cbind(cell_line_CNA_high_level_only[genes_in_all_4_files, cell_lines_with_both_MUT_and_CNA], tumor_CNA_high_level_only[genes_in_all_4_files, tumors_with_both_MUT_and_CNA])
  composite_MUT <- cbind(cell_line_MUT[genes_with_MUT_in_both, cell_lines_with_both_MUT_and_CNA], tumor_MUT[genes_with_MUT_in_both, tumors_with_both_MUT_and_CNA])
  
  rownames(composite_MUT) <- paste(rownames(composite_MUT), "MUT", sep = "_")
  rownames(composite_CNA) <- paste(rownames(composite_CNA), "CNA", sep = "_")
  rownames(composite_CNA_high_level_only) <- paste(rownames(composite_CNA_high_level_only), "CNA", sep = "_")
  
  # Generate matrix and convert to matrix 
  composite_mat <- rbind(composite_MUT, composite_CNA)
  composite_mat <- as.matrix(composite_mat)
  
  # WRITE COMPOSITE
  write.table(composite_mat,
              file = output_composite_alteration_matrix_file,
              sep = "\t",
              quote = FALSE)
  #composite_mat_high_level_only <- rbind(composite_MUT,composite_CNA_high_level_only)
  
  # Calculation of alteration frequencies
  # Assign frequency weights as (freq. of alteration of gene)/(mean freq. of alteration across all genes) - "rewarding recurrent changes"
  overall_alt_freq <- length(which((composite_mat[]) != 0)) / (length(which((composite_mat[]) == 0)) + length(which((composite_mat[]) !=
                                                                                        0)))
  freq_alt <- rep(0, nrow(composite_mat))
  freq_alt <- apply(composite_mat, 1, compute_freq_alt)
  #freq_alt_high_level  <- apply(composite_mat_high_level_only,1,compute_freq_alt)
  freq_alt_mut_tumors <- apply(composite_MUT[, tumors_with_both_MUT_and_CNA], 1, compute_freq_alt)
  freq_alt_cna_tumors <- apply(composite_CNA[, tumors_with_both_MUT_and_CNA], 1, compute_freq_alt)
  
  freq_alt_samplewise <- apply(composite_mat, 2, compute_freq_alt)
  #freq_alt_samplewise_CNA_high_level_only <- apply(composite_mat_high_level_only,2,compute_freq_alt)
  
  composite_mat <- composite_mat[, which(freq_alt_samplewise > 0)]
  #composite_mat_high_level_only <- composite_mat_high_level_only[,which(freq_alt_samplewise > 0)]
  
  names(freq_alt) <- rownames(composite_mat)
  freq_weights <- rep(1, nrow(composite_mat))
  freq_weights <- freq_alt / overall_alt_freq
  names(freq_weights) <- rownames(composite_mat)
  
  # GET WEIGHTS ----
  # Read in user-provided weights
  known_cancer_genes_and_weights_all <-
    read.table(
      pancancer_gene_weights_file,
      sep = "\t",
      header = TRUE,
      row.names = 1
    )
  known_cancer_genes_and_weights <- as.matrix(known_cancer_genes_and_weights_all[intersect(rownames(known_cancer_genes_and_weights_all), rownames(composite_mat)), ]) # To eliminate entries not present in alteration matrix, if any
  genes_and_weights_all <-
    read.table(
      cancer_specific_gene_weights_file,
      sep = "\t",
      header = TRUE,
      row.names = 1
    )
  # To eliminate entries not present in alteration matrix, if any
  genes_and_weights <- as.matrix(genes_and_weights_all[intersect(rownames(genes_and_weights_all), rownames(composite_mat)), ]) 
  rownames(genes_and_weights) <- intersect(rownames(genes_and_weights_all), rownames(composite_mat))
  rownames(known_cancer_genes_and_weights) <- intersect(rownames(known_cancer_genes_and_weights_all), rownames(composite_mat))
  
  # FIX WHAT IS THIS ? 
  annotation_weights <- rep(NA, nrow(composite_mat))
  annotation_weights[grep("_CNA", rownames(composite_mat))] <- CNA_default_weight
  annotation_weights[grep("_MUT", rownames(composite_mat))] <- MUT_default_weight
  
  names(annotation_weights) <- rownames(composite_mat)
  for (i in 1:nrow(known_cancer_genes_and_weights))
    # Overwrite default weights for known cancer genes
    annotation_weights[rownames(known_cancer_genes_and_weights)[i]] = known_cancer_genes_and_weights[i, ]
  for (i in 1:nrow(genes_and_weights))
    # Overwrite default weights for input provided weights
    annotation_weights[rownames(genes_and_weights)[i]] = genes_and_weights[i, ]
  
  gene_weights <- rep(1, nrow(composite_mat))
  names(gene_weights) <- rownames(composite_mat)
  
  gene_weights <- annotation_weights # if using user-provided weights only
  gene_weights <- gene_weights / max(gene_weights) # map to 0-1
  
  cor_unweighted <- cor(composite_mat)
  
  switch(distance_similarity_measure, 
         weighted_correlation={
           # CALCULATE CORRELATIONS ----
           # Including low-level CNAs
           cor_weighted <- calc_weighted_corr(as.matrix(composite_mat),
                                              as.matrix(composite_mat),
                                              gene_weights)
           # Excluding low-levels CNAs
           #cor_weighted_high_level_only <- calc_weighted_corr(as.matrix(composite_mat_high_level_only),as.matrix(composite_mat_high_level_only),gene_weights)
           
          # Convert to distance, and call multidimensional scaling via isoMDS
          dist <- 1 - as.matrix(cor_weighted)
          isomdsfit <- isoMDS(dist, k=2)
         },
         generalized_jaccard={
           # Calculate weighted distance based on Jaccard's coefficient
           weighted_distance_excluding_zero_zero_matches <- apply(composite_mat, 2, function(x_i,weights=gene_weights) 
             sapply(1:ncol(composite_mat), function(j) pair_dist(x_i, composite_mat[,j], gene_weights))) # repeatedly apply function for weighted distance between a pair of coulmns/vectors
           
           # Change missing or small values
           weighted_distance_excluding_zero_zero_matches[which(is.na(weighted_distance_excluding_zero_zero_matches))] <- 0
           weighted_distance_excluding_zero_zero_matches <- weighted_distance_excluding_zero_zero_matches + 1e-6
           # Call multidimensional scaling via isoMDS
           dist <- weighted_distance_excluding_zero_zero_matches
           isomdsfit <-  isoMDS(dist, k=2)
         }
  )

  results <- list(
    dist = dist,
    isomdsfit = isomdsfit, 
    cor_unweighted = cor_unweighted,
    composite_mat = composite_mat,
    cell_lines_with_both_MUT_and_CNA = cell_lines_with_both_MUT_and_CNA,
    tumors_with_both_MUT_and_CNA = tumors_with_both_MUT_and_CNA
  )
  
  return(results)
}