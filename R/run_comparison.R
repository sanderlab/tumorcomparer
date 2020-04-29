#' Run a comparison between between two cohorts (e.g. cell lines and tumors)
#'
#' @param available_data_types a vector of data types to be analyzed; the order of these much match distance_similarity_measures
#' 
#' @param cna_data_type_weight TODO
#' @param mut_data_type_weight TODO
#' @param exp_data_type_weight TODO
#' 
#' @param cna_default_weight default (background) weight for copy number alterations (CNA) (DEFAULT: 0.01). 
#'   Default weights are assigned to genes not known to be important in the specific cancer type or cancer in general
#' @param mut_default_weight default (background) weight for mutation alterations (MUT) (DEFAULT: 0.01); CNA_default_weight
#' @param exp_default_weight default (background) weight for mRNA gene expression values (EXP) (DEFAULT: 0.01)
#' 
#' @param cna_known_cancer_gene_weight a default weight for genes important in cancer (DEFAULT: 0.1); 
#'   an example would genes in the TCGA pan-cancer CNA list (i.e. GISTIC peaks). Known cancer gene weights are 
#'   default weight for genes important in cancer (an example would genes derived from TCGA pan-cancer analyses or literature)
#' @param mut_known_cancer_gene_weight a default weight for genes important in cancer (DEFAULT: 0.1);
#'   an example would genes in the TCGA pan-cancer mutation list (i.e. MUTSIG SMGs)
#' @param exp_known_cancer_gene_weight a default weight for genes important in cancer(DEFAULT: 0.1); 
#'   an example would genes in the COSMIC Cancer Gene Census
#' 
#' @param tumor_mut_file a file with binary mutation data for tumors 
#' @param tumor_cna_file a file with 5-valued GISTIC data for tumors
#' @param tumor_exp_file a file with gene expression data for tumors
#' 
#' @param cell_line_mut_file a file with binary mutation data for cell lines
#' @param cell_line_cna_file a file with 5-valued GISTIC data for cell lines
#' @param cell_line_exp_file a file with gene expression data for cell lines
#' 
#' @param known_cancer_gene_weights_mut_file a file with weights for genes known
#'   to be recurrently altered/mutated in cancer (e.g. recurrently mutated genes in TCGA pan-cancer analyses). 
#'   A two-column tab-delimited file - the first column has the gene names and the second column specifies the weights.
#' @param known_cancer_gene_weights_cna_file see known_cancer_gene_weights_mut_file
#' @param known_cancer_gene_weights_exp_file see known_cancer_gene_weights_mut_file
#' 
#' @param cancer_specific_gene_weights_mut_file a file with weights for cancer-specific
#'   set of recurrently mutated genes. A tab-delimited file - the first column has the gene names,
#'   and the second column specifies the weights.
#' @param cancer_specific_gene_weights_cna_file see cancer_specific_gene_weights_mut_file 
#' @param cancer_specific_gene_weights_exp_file see cancer_specific_gene_weights_mut_file 
#' 
#' @param distance_similarity_measures a named vector of distance/similarity measures: see details. 
#'   OPTIONS: weighted_correlation and generalized_jaccard - must be in the order mut, cna, exp. Currently, 
#'   generalized_jaccard is used for mut and cna data, and weighted_correlation for exp data
#'
#' @return a list with multiple items. NOTE: The values of the dist and isomdsfit will depend on parameter, distance_similarity_measure.
#' * dist_mat: a matrix of pairwise distances
#' * isomdsfit: a two-column (2-dimension) fitting of the distances reduced to two dimensions via MDS - multidimensional scaling, using the isoMDS function
#' * cor_unweighted: a matrix of unweighted pairwise correlations
#' * composite_mat: the composite matrix (see Details)
#' * cell_line_ids: a vector of cell line IDs/names with all data types
#' * tumor_ids: a vector of tumor IDs with all data types
#' @md
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom MASS isoMDS
run_comparison <- function(available_data_types=c("mut", "cna", "exp"), 
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
                           distance_similarity_measures=c("generalized_jaccard", "weighted_correlation", "weighted_correlation")
                           #distance_similarity_measures=c("generalized_jaccard", "generalized_jaccard", "weighted_correlation")
                           ) {
  
  # CHECK available_data_types and distance_similarity_measures ----
  if(length(available_data_types) < 1) {
    stop("ERROR: At least one data type: mut, cna, exp must be provided for available_data_types")  
  } 
  
  if(length(available_data_types) != length(distance_similarity_measures)) {
    stop("ERROR: length(available_data_types) must equal length(distance_similarity_measures)")  
  } 
       
  # LOAD DATA ----            
  isomdsfit_by_data_type <- list() 
  dist_mat_by_data_type <- list() 
  composite_mat_by_data_type <- list() 
  gene_weights_by_data_type <- list() 
  cancer_specific_gene_weights_by_data_type <- list()
  known_cancer_gene_weights_by_data_type <- list()
  count <- 1
  
  for(data_type in available_data_types) {
    #cat("DEBUG: ", data_type, "\n")
    
    if(data_type == "mut") {
      #cat("DEBUG\n") 
      
      mut <- generate_composite_mat_and_gene_weights(
        default_weight=mut_default_weight,
        known_cancer_gene_weight=mut_known_cancer_gene_weight,
        tumor_file=tumor_mut_file,
        cell_line_file=cell_line_mut_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_mut_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_mut_file,
        is_discrete=TRUE,
        distance_similarity_measure=distance_similarity_measures[count])
      
      isomdsfit_by_data_type[["mut"]] <- mut$isomdsfit
      dist_mat_by_data_type[["mut"]] <- mut$dist_mat
      composite_mat_by_data_type[["mut"]] <- mut$composite_mat
      gene_weights_by_data_type[["mut"]] <- mut$gene_weights
      known_cancer_gene_weights_by_data_type[["mut"]] <- mut$known_cancer_gene_weights
      cancer_specific_gene_weights_by_data_type[["mut"]] <- mut$cancer_specific_gene_weights
    }
    
    if(data_type == "cna") {
      cna <- generate_composite_mat_and_gene_weights(
        default_weight=cna_default_weight,
        known_cancer_gene_weight=cna_known_cancer_gene_weight,
        tumor_file=tumor_cna_file,
        cell_line_file=cell_line_cna_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_cna_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_cna_file, 
        is_discrete=FALSE,
        distance_similarity_measure=distance_similarity_measures[count])
      
      isomdsfit_by_data_type[["cna"]] <- cna$isomdsfit
      dist_mat_by_data_type[["cna"]] <- cna$dist_mat
      composite_mat_by_data_type[["cna"]] <- cna$composite_mat
      gene_weights_by_data_type[["cna"]] <- cna$gene_weights
      known_cancer_gene_weights_by_data_type[["cna"]] <- cna$known_cancer_gene_weights
      cancer_specific_gene_weights_by_data_type[["cna"]] <- cna$cancer_specific_gene_weights
    }
    
    if(data_type == "exp") {
      exp <- generate_composite_mat_and_gene_weights(
        default_weight=exp_default_weight,
        known_cancer_gene_weight=exp_known_cancer_gene_weight,
        tumor_file=tumor_exp_file,
        cell_line_file=cell_line_exp_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_exp_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_exp_file,
        is_discrete=FALSE,
        distance_similarity_measure=distance_similarity_measures[count])
      
      isomdsfit_by_data_type[["exp"]] <- exp$isomdsfit
      dist_mat_by_data_type[["exp"]] <- exp$dist_mat
      composite_mat_by_data_type[["exp"]] <- exp$composite_mat
      gene_weights_by_data_type[["exp"]] <- exp$gene_weights
      known_cancer_gene_weights_by_data_type[["exp"]] <- exp$known_cancer_gene_weights
      cancer_specific_gene_weights_by_data_type[["exp"]] <- exp$cancer_specific_gene_weights
    }
    
    count <- count + 1
  }
  
  # SUM DATA WEIGHTS ----
  sum_data_type_weights <- 0 
  
  for(data_type in available_data_types) {
    if(data_type == "mut") {
      sum_data_type_weights <- sum_data_type_weights + mut_data_type_weight
    }
    
    if(data_type == "cna") {
      sum_data_type_weights <- sum_data_type_weights + cna_data_type_weight
    }
    
    if(data_type == "exp") {
      sum_data_type_weights <- sum_data_type_weights + exp_data_type_weight
    }
  }

  # Check for weights sum
  if(sum_data_type_weights != 1) {
    stop("ERROR: Sum of *_data_type_weights must sum up to 1")  
  } 
  
  # INITIALIZE ALL NECESSARY LISTS ----
  combined_samples_list <- list() 
  combined_tumor_ids_list <- list() 
  combined_cell_line_ids_list <- list() 
  
  for(data_type in available_data_types) {
    if(data_type == "mut") {
      mut_samples <- colnames(mut$dist_mat) 
      mut_tumor_ids <- mut$tumor_ids 
      mut_cell_line_ids <- mut$cell_line_ids
      
      combined_samples_list[["mut"]] <- mut_samples
      combined_tumor_ids_list[["mut"]] <- mut_tumor_ids
      combined_cell_line_ids_list[["mut"]] <- mut_cell_line_ids
    }
    
    if(data_type == "cna") {
      cna_samples <- colnames(cna$dist_mat) 
      cna_tumor_ids <- cna$tumor_ids
      cna_cell_line_ids <- cna$cell_line_ids
      
      combined_samples_list[["cna"]] <- cna_samples
      combined_tumor_ids_list[["cna"]] <- cna_tumor_ids
      combined_cell_line_ids_list[["cna"]] <- cna_cell_line_ids
    }
    
    if(data_type == "exp") {
      exp_samples <- colnames(exp$dist_mat) 
      exp_tumor_ids <- exp$tumor_ids
      exp_cell_line_ids <- exp$cell_line_ids
      
      combined_samples_list[["exp"]] <- exp_samples
      combined_tumor_ids_list[["exp"]] <- exp_tumor_ids
      combined_cell_line_ids_list[["exp"]] <- exp_cell_line_ids
    }
  }
  
  combined_samples <- Reduce(intersect, combined_samples_list)
  combined_tumor_ids <- Reduce(intersect, combined_tumor_ids_list)
  combined_cell_line_ids <- Reduce(intersect, combined_cell_line_ids_list)
  
  # CALCULATE COMBINED_DIST AND ISOMDS ----
  # Idea: If null, set to the matrix check, otherwise add to the existing combined_dist
  combined_dist <- NULL 
  
  for(data_type in available_data_types) {
    if(data_type == "mut") {
      if(is.null(combined_dist)) {
        combined_dist <- mut_data_type_weight*mut$dist_mat[combined_samples, combined_samples]
      } else {
        combined_dist <- combined_dist + mut_data_type_weight*mut$dist_mat[combined_samples, combined_samples]
      }
    }
    
    if(data_type == "cna") {
      if(is.null(combined_dist)) {
        combined_dist <- cna_data_type_weight*cna$dist_mat[combined_samples, combined_samples]
      } else {
        combined_dist <- combined_dist + cna_data_type_weight*cna$dist_mat[combined_samples, combined_samples]
      }
    }
    
    if(data_type == "exp") {
      if(is.null(combined_dist)) {
        combined_dist <- exp_data_type_weight*exp$dist_mat[combined_samples, combined_samples]
      } else {
        combined_dist <- combined_dist + exp_data_type_weight*exp$dist_mat[combined_samples, combined_samples]
      }
    }
  }
  
  a <- 1
  
  # RUN ISOMDS ----
  isomdsfit <-  isoMDS(combined_dist, k=2)  
  
  # MERGE RESULTS ----
  results <- list(
    dist_mat = combined_dist,
    dist_mat_by_data_type = dist_mat_by_data_type,
    composite_mat_by_data_type = composite_mat_by_data_type,
    gene_weights_by_data_type = gene_weights_by_data_type,
    isomdsfit = isomdsfit,
    isomdsfit_by_data_type = isomdsfit_by_data_type,
    cell_line_ids = combined_cell_line_ids,
    tumor_ids = combined_tumor_ids,
    known_cancer_gene_weights_by_data_type = known_cancer_gene_weights_by_data_type,
    cancer_specific_gene_weights_by_data_type = cancer_specific_gene_weights_by_data_type
  )
  
  return(results)
}

