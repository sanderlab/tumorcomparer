#' Run a comparison between between two cohorts (e.g. cell lines and tumors)
#'
#' @param data_types a vector of data types to be analyzed; the order of these much match distance_similarity_measures
#' 
#' @param CNA_default_weight default (background) weight for copy number alterations (CNA) (DEFAULT: 0.01). 
#'   Default weights are assigned to genes not known to be important in the specific cancer type or cancer in general
#' @param MUT_default_weight default (background) weight for mutation alterations (MUT) (DEFAULT: 0.01); CNA_default_weight
#' @param EXP_default_weight default (background) weight for mRNA gene expression values (EXP) (DEFAULT: 0.01)
#' 
#' @param CNA_known_cancer_gene_weight a default weight for genes important in cancer (DEFAULT: 0.1); 
#'   an example would genes in the TCGA pan-cancer CNA list (i.e. GISTIC peaks). Known cancer gene weights are 
#'   default weight for genes important in cancer (an example would genes derived from TCGA pan-cancer analyses or literature)
#' @param MUT_known_cancer_gene_weight a default weight for genes important in cancer (DEFAULT: 0.1);
#'   an example would genes in the TCGA pan-cancer mutation list (i.e. MUTSIG SMGs)
#' @param EXP_known_cancer_gene_weight a default weight for genes important in cancer(DEFAULT: 0.1); 
#'   an example would genes in the COSMIC Cancer Gene Census
#' 
#' @param tumor_mut_file a file with binary mutation data for tumors 
#' @param tumor_cna_file a file with 5-valued GISTIC data for tumors
#' @param tumor_exp_file a file with gene expression data for tumors
#' @param cell_line_mut_file a file with binary mutation data for cell lines
#' @param cell_line_cna_file a file with 5-valued GISTIC data for cell lines
#' @param cell_line_exp_file a file with gene expression data for cell lines

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
#' @param distance_similarity_measures a named vector of distance/similarity measures (See Details) 
#'   (OPTIONS: "weighted_correlation", "generalized_jaccard" - must be in the order mut, cna, exp) 
#'   currently, "generalized_jaccard" is used for mut and cna data, and "weighted_correlation" for exp data
#' 
#' 
#'    
#' @return a list with multiple items. NOTE: The values of the dist and isomdsfit will
#'  depend on parameter "distance_similarity_measure".
#' \itemize{
#'   \item{"dist_mat"}{a matrix of pairwise distances}
#'   \item{"isomdsfit"}{a two-column (2-dimension) fitting of the distances reduced to 
#'   two dimensions via MDS - multidimensional scaling, using the isoMDS function}
#'   \item{"cor_unweighted"}{a matrix of unweighted pairwise correlations}
#'   \item{"composite_mat"}{the composite matrix (see Details)}
#'   \item{"cell_lines_with_both_MUT_and_CNA"}{a vector of cell line IDs/names with both mutation
#'    (MUT) and copy number alteration (CNA) information}
#'   \item{"tumors_with_both_MUT_and_CNA"}{a vector of tumor IDs with both mutation 
#'   (MUT) and copy number alteration (CNA) information}
#' }
#'
#' @author Rileen Sinha (rileen@gmail.com), Augustin Luna (aluna@jimmy.harvard.edu)
#'
#' @concept tumorcomparer
#' @export
#' 
#' @importFrom MASS isoMDS
#' @importFrom cluster daisy
#' @importFrom utils read.table write.table
#' @importFrom stats cor
run_comparison <- function(data_types=c("mut", "cna", "exp"), 
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
                           cancer_specific_gene_weights_exp_file="Genes_and_weights_exp.txt"), 
                           distance_similarity_measures=c("generalized jaccard", "generalized jaccard", "weighted correlation")
                           ) {
                           
  for(data_type in available_data_types) {
    if(names(data_type) == "mut") {
      mut <- generate_composite_mat_and_gene_weights(
        default_weight=mut_default_weight,
        known_cancer_gene_weight=mut_known_cancer_gene_weight,
        tumor_file=tumor_mut_file,
        cell_line_file=cell_line_mut_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_mut_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_mut_file,
        is_discrete=TRUE,
        distance_similarity_measure=data_type)
    }
    
    if(names(data_type) == "cna") {
      cna <- generate_composite_mat_and_gene_weights(
        default_weight=cna_default_weight,
        known_cancer_gene_weight=cna_known_cancer_gene_weight,
        tumor_file=tumor_cna_file,
        cell_line_file=cell_line_cna_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_cna_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_cna_file, 
        is_discrete=TRUE,
        distance_similarity_measure=data_type)
    }
    
    if(names(data_type) == "exp") {
      exp <- generate_composite_mat_and_gene_weights(
        default_weight=exp_default_weight,
        known_cancer_gene_weight=exp_known_cancer_gene_weight,
        tumor_file=tumor_exp_file,
        cell_line_file=cell_line_exp_file,
        known_cancer_gene_weights_file=known_cancer_gene_weights_exp_file,
        cancer_specific_gene_weights_file=cancer_specific_gene_weights_exp_file,
        is_discrete=FALSE,
        distance_similarity_measure=data_type)
    }
  }
  
  # Check data type weights
  sum_data_type_weights <- mut_data_type_weight + cna_data_type_weight + exp_data_type_weight

  if(sum_data_type_weights != 1) {
    stop("ERROR: Sum of *_data_type_weights must sum up to 1")  
  } 
  
  mut_samples <- colnames(mut$dist_mat) 
  cna_samples <- colnames(cna$dist_mat) 
  exp_samples <- colnames(exp$dist_mat) 
  combined_samples <- Reduce(intersect, list(mut_samples, cna_samples, exp_samples))
  
  combined_dist <- mut_data_type_weight*mut$dist_mat[combined_samples, combined_samples] + 
    cna_data_type_weight*cna$dist_mat[combined_samples, combined_samples] + 
    exp_data_type_weight*exp$dist_mat[combined_samples, combined_samples]

  # Calculate isomdsfit 
  isomdsfit <-  isoMDS(combined_dist, k=2)  
  
  mut_tumor_ids <- colnames(mut$tumor_ids) 
  cna_tumor_ids <- colnames(cna$tumor_ids) 
  exp_tumor_ids <- colnames(exp$tumor_ids) 
  combined_tumor_ids <- Reduce(intersect, list(mut_tumor_ids, cna_tumor_ids, exp_tumor_ids))
  
  mut_cell_line_ids <- colnames(mut$cell_line_ids) 
  cna_cell_line_ids <- colnames(cna$cell_line_ids) 
  exp_cell_line_ids <- colnames(exp$cell_line_ids) 
  combined_cell_lines_ids <- Reduce(intersect, list(mut_cell_line_ids, cna_cell_line_ids, exp_cell_line_ids))
  
  results <- list(
    dist_mat = combined_dist,
    isomdsfit = isomdsfit, 
    isomdsfit_by_data_type=list(mut=mut$isomdsfit, cna=cna$isomdsfit, exp=exp$isomdsfit)
    cell_line_ids = combined_cell_line_ids,
    tumor_ids = combined_tumor_ids
  )
  
  return(results)
}

