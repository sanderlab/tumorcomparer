#' Run a comparison between between two cohorts (e.g. cell lines and tumors)
#' 
#' @param default_weight see run_comparison
#' @param tumor_file see run_comparison
#' @param cell_line_file see run_comparison
#' @param known_cancer_gene_weights_file see run_comparison
#' @param cancer_specific_gene_weights_file see run_comparison
#' 
#' @details The composite matrix is a single matrix where the columns are samples 
#'   (i.e. tumors AND cell line IDs) and the rows are an rbind() of mutations 
#'   (with 1 or 0 outputs for each sample), copy number alterations from 
#'   GISTIC (with values -2, -1, 0, 1, 2), or gene expression values. Available similarity/distance measures include: 
#'   \itemize{
#'     \item{"weighted_correlation"}{Weighted correlation, based on weighted means and standard deviations}
#'     \item{"generalized_jaccard"}{A weighted distance based on the Jaccard coefficient}
#'  }
#' 
#' @return a list with multiple items. NOTE: The values of the dist and isomdsfit will
#'  depend on parameter "distance_similarity_measure".
#' \itemize{
#'   \item{"dist_mat"}{a matrix of pairwise distances}
#'   \item{"isomdsfit"}{a two-column (2-dimension) fitting of the distances reduced to 
#'   two dimensions via MDS - multidimensional scaling, using the isoMDS function}
#'   \item{"cor_unweighted"}{a matrix of unweighted pairwise correlations}
#'   \item{"composite_mat"}{the composite matrix (see Details)}
#'   \item{"cell_lines_ids"}{a vector of cell line IDs/names}
#'   \item{"tumors_ids"}{a vector of tumor IDs}
#' }
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#' 
#' @examples 
#' mut_default_weight <- 0.01
#' tumor_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "tumor_mut.txt", 
#'   package="tumorcomparer")
#' cell_line_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "cell_line_mut.txt", 
#'   package="tumorcomparer")
#' known_cancer_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", 
#'   "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
#' cancer_specific_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", 
#'   "Genes_and_weights_mut.txt", package="tumorcomparer")
#' 
#' mut <- generate_composite_mat_and_gene_weights(
#'   default_weight=mut_default_weight,
#'   tumor_file=tumor_mut_file,
#'   cell_line_file=cell_line_mut_file,
#'   known_cancer_gene_weights_file=known_cancer_gene_weights_mut_file,
#'   cancer_specific_gene_weights_file=cancer_specific_gene_weights_mut_file)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom MASS isoMDS
#' @importFrom cluster daisy
#' @importFrom utils read.table write.table
#' @importFrom stats cor
generate_composite_mat_and_gene_weights <- function(default_weight, 
                                                    tumor_file, 
                                                    cell_line_file, 
                                                    known_cancer_gene_weights_file, 
                                                    cancer_specific_gene_weights_file, 
                                                    distance_similarity_measure) {

  # GET INTERSECTING GENES BETWEEN TUMORS AND CELL LINES ----
  tumor <- read.table(tumor_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  cell_line <- read.table(cell_line_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  tumor_ids <- colnames(tumor)
  cell_line_ids <- colnames(cell_line) 
  
  # Genes in both tumors and cell
  genes_in_both <- intersect(rownames(tumor), rownames(cell_line))

  composite_mat <- cbind(cell_line[genes_in_both,], tumor[genes_in_both,])

  # Check that tumor and cell_line data are either both discrete or continuous 
  # and set the distance_similarity_measure to run
  if(all((composite_mat %% 1) == 0)) {
    is_discrete <- TRUE
    distance_similarity_measure <- "generalized_jaccard"      
  } else {
    is_discrete <- FALSE 
    distance_similarity_measure <- "weighted_correlation"      
  }

  # For discrete data (mut, cna) - compute alteration frequencies, remove samples which have no alterations
  if(is_discrete) { 
    # Calculation of alteration frequencies
    # Assign frequency weights as (freq. of alteration of gene)/(mean freq. of alteration across all genes) - "rewarding recurrent changes"
    overall_alt_freq <- length(which((composite_mat[]) != 0)) / (length(which((composite_mat[]) == 0)) + length(which((composite_mat[]) != 0)))
  
    freq_alt <- rep(0, nrow(composite_mat))
    freq_alt <- apply(composite_mat, 1, compute_freq_alt)
    freq_alt_tumors <- apply(composite_mat[, colnames(tumor)], 1, compute_freq_alt)
    freq_alt_cell_lines <- apply(composite_mat[, colnames(cell_line)], 1, compute_freq_alt)
    
    freq_alt_samplewise <- apply(composite_mat, 2, compute_freq_alt)
  
    #composite_mat <- composite_mat[, which(freq_alt_samplewise > 0)]
    
    names(freq_alt) <- rownames(composite_mat)
    
    freq_weights <- rep(1, nrow(composite_mat))
    freq_weights <- freq_alt / overall_alt_freq
    names(freq_weights) <- rownames(composite_mat)
  }
 
  #Intitalize weights with default weights
  gene_weights <- rep(default_weight, nrow(composite_mat))
  names(gene_weights) <- rownames(composite_mat)

 
  # GET WEIGHTS ----
  # Read in user-provided weights for known cancer genes
  known_cancer_genes_and_weights_all <-
    read.table(
      known_cancer_gene_weights_file,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE
    )
  rownames(known_cancer_genes_and_weights_all) <- trimws(rownames(known_cancer_genes_and_weights_all)) # trim whitespace, if any

  gene_weights[intersect(names(gene_weights), rownames(known_cancer_genes_and_weights_all))] <- known_cancer_genes_and_weights_all[intersect(names(gene_weights),rownames(known_cancer_genes_and_weights_all)),1]


  # Read in user-provided weights for cancer-specific genes
  genes_and_weights_all <-
    read.table(
      cancer_specific_gene_weights_file,
      sep = "\t",
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE
    )
  rownames(genes_and_weights_all) <- trimws(rownames(genes_and_weights_all)) # trim whitespace, if any
    
  gene_weights[intersect(names(gene_weights), rownames(genes_and_weights_all))] <- genes_and_weights_all[intersect(names(gene_weights),rownames(genes_and_weights_all)),1]
  
  gene_weights <- gene_weights / max(gene_weights + 1e-6) # map to 0-1
  
  #cor_unweighted <- cor(composite_mat)  
  
  if(distance_similarity_measure == "weighted_correlation") {
    # CALCULATE CORRELATIONS ----
    # Including low-level CNAs
    cor_weighted <- calc_weighted_corr(as.matrix(composite_mat),
                                       as.matrix(composite_mat),
                                       gene_weights)
    cor_weighted[which(is.na(cor_weighted))] <- 0
    cor_weighted <- cor_weighted - 1e-6
    # Excluding low-levels CNAs
    #cor_weighted_high_level_only <- calc_weighted_corr(as.matrix(composite_mat_high_level_only),as.matrix(composite_mat_high_level_only),gene_weights)
    
    # Convert to distance, and call multidimensional scaling via isoMDS
    dist_mat <- 1 - as.matrix(cor_weighted)
    isomdsfit <- isoMDS(dist_mat, k=2)
  } else if(distance_similarity_measure == "generalized_jaccard") {
    # Calculate weighted distance based on Jaccard's coefficient
    weighted_distance_excluding_zero_zero_matches <- apply(composite_mat, 2, function(x_i,weights=gene_weights) 
      sapply(1:ncol(composite_mat), function(j) pair_dist(x_i, composite_mat[,j], gene_weights))) # repeatedly apply function for weighted distance between a pair of coulmns/vectors
    
    # Change missing or small values
    weighted_distance_excluding_zero_zero_matches[which(is.na(weighted_distance_excluding_zero_zero_matches))] <- 0
    weighted_distance_excluding_zero_zero_matches <- weighted_distance_excluding_zero_zero_matches + 1e-6
    # Call multidimensional scaling via isoMDS
    dist_mat <- weighted_distance_excluding_zero_zero_matches
    rownames(dist_mat) <- colnames(composite_mat)
    colnames(dist_mat) <- colnames(composite_mat)
    isomdsfit <-  isoMDS(dist_mat, k=2)  
  } else {
    stop("ERROR: Unknown distance_similarity_measure: ", distance_similarity_measure)
  }
 
  results <- list(
    dist_mat = dist_mat,
    isomdsfit = isomdsfit, 
    composite_mat = composite_mat,
    gene_weights = gene_weights,
    cancer_specific_gene_weights = genes_and_weights_all,
    known_cancer_gene_weights = known_cancer_genes_and_weights_all,
    cell_line_ids = cell_line_ids,
    tumor_ids = tumor_ids
  )
  
  return(results)
}

