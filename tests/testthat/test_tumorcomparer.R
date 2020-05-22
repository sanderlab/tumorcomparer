# All tests are done on files in package using system.file()

test_that("sanity_check", {
  expect_equal(2 * 2, 4)
})

test_that("return_first_part", {
  tmp <- return_first_part("22RV1_PROSTRATE")
  expect_identical(tmp, "22RV1")
})

test_that("keep_only_high_level_cnas", {
  set.seed(1)
  tmp <- sample(c(-2, -1, 0, 1, 2), 100, replace=TRUE)
  mat <- matrix(tmp, 10, 10)
  
  results <- as.vector(keep_only_high_level_cnas(cna_mat=mat))
  expect_false(any(results == 1 | results == -1))
})

test_that("keep_only_high_level_cnas", {
  t1 <- c(1,1,1,1,1,1,2,2,2,3,3,3,4,4)
  t2 <- rev(t1)
  weight <- c(.5,.5,.5,.5,.5,1,1,1,1,2,2,2,2,2)
  
  results <- calc_weighted_corr(t1, t2, weight)

  expect_equal(results, matrix(-0.8108894, 1, 1))
})

test_that("compute_freq_alt", {
  set.seed(1)
  tmp <- sample(c(0, 1), 100, replace=TRUE)
  mat <- matrix(tmp, 10, 10)
  expect_equal(0.51, compute_freq_alt(mat))
})

test_that("categorize_cell_lines", {
  comparison_result <- readRDS(system.file("test_output", "ov_comparison_result.rds", package="tumorcomparer"))
  
  categorization_list <- categorize_cell_lines(
    num_tumors_for_comparison=length(comparison_result$tumor_ids)-1, 
    dist_mat=comparison_result$dist_mat,
    cell_line_ids=comparison_result$cell_line_ids,
    tumor_ids=comparison_result$tumor_ids,
    trim_cell_line_names=FALSE) 
  
  expect_equal(categorization_list$categorization$Sample_ID[c(1,6,13)], c("OAW-28", "OV-17R", "OVK-18"))
  expect_equal(categorization_list$categorization$Category[c(1,6,13)], c("Moderately Good", "Poor", "Outlier"))
})

test_that("generate_composite_mat_and_gene_weights", {
  expect_equal(TRUE, FALSE)
})

test_that("map_mean_similarity_to_gradient", {
  comparison_result <- readRDS(system.file("test_output", "ov_comparison_result.rds", package="tumorcomparer"))
  
  categorization_list <- categorize_cell_lines(
    num_tumors_for_comparison=length(comparison_result$tumor_ids)-1, 
    dist_mat=comparison_result$dist_mat,
    cell_line_ids=comparison_result$cell_line_ids,
    tumor_ids=comparison_result$tumor_ids,
    trim_cell_line_names=FALSE) 
  
  result <- map_mean_similarity_to_gradient(
    mean_similarity_cell_line_to_k_nearest_tumors=categorization_list$mean_similarity_cell_line_to_k_nearest_tumors,
    mean_similarity_tumor_to_k_nearest_tumors=categorization_list$mean_similarity_tumor_to_k_nearest_tumors,
    col1="orange",
    col2="blue", 
    numshades=100)
  
  expect_equal(result, c("#DD8F21", "#F9A105", "#C37E3B"))
})

test_that("pair_dist", {
  set.seed(1)
  
  n <- 100
  x <- sample(c(0, 1), n, replace=TRUE)
  y <- sample(c(0, 1), n, replace=TRUE)
  weights <- rnorm(n)
  result <- pair_dist(x, y, weights)
  
  expect_equal(result, 2.398, tolerance = 0.01)
})

test_that("plot_mds", {
  expect_equal(TRUE, FALSE)
})

test_that("run_comparison", {
  tumor_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "tumor_mut.txt", package="tumorcomparer")
  tumor_cna_file <- system.file("extdata", "ovarian_tcga_cclp", "tumor_cna.txt", package="tumorcomparer")
  tumor_exp_file <- system.file("extdata", "ovarian_tcga_cclp", "tumor_exp.txt", package="tumorcomparer")
  
  cell_line_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "cell_line_mut.txt", package="tumorcomparer")
  cell_line_cna_file <- system.file("extdata", "ovarian_tcga_cclp", "cell_line_cna.txt", package="tumorcomparer")
  cell_line_exp_file <- system.file("extdata", "ovarian_tcga_cclp", "cell_line_exp.txt", package="tumorcomparer")
  
  known_cancer_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
  known_cancer_gene_weights_cna_file <- system.file("extdata", "ovarian_tcga_cclp", "default_weights_for_known_cancer_genes_cna.txt", package="tumorcomparer")
  known_cancer_gene_weights_exp_file <- system.file("extdata", "ovarian_tcga_cclp", "default_weights_for_known_cancer_genes_exp.txt", package="tumorcomparer")
  
  cancer_specific_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "Genes_and_weights_mut.txt", package="tumorcomparer")
  cancer_specific_gene_weights_cna_file <- system.file("extdata", "ovarian_tcga_cclp", "Genes_and_weights_cna.txt", package="tumorcomparer")
  cancer_specific_gene_weights_exp_file <- system.file("extdata", "ovarian_tcga_cclp", "Genes_and_weights_exp.txt", package="tumorcomparer")
  
  comparison_result <- run_comparison(
    available_data_types=c("mut", "cna", "exp"), 
    mut_data_type_weight = 1/3,
    cna_data_type_weight = 1/3,
    exp_data_type_weight = 1/3,
    cna_default_weight=0.01, 
    mut_default_weight=0.01,
    exp_default_weight=0.01,
    cna_known_cancer_gene_weight=0.1, 
    mut_known_cancer_gene_weight=0.1, 
    exp_known_cancer_gene_weight=0.1, 
    tumor_mut_file=tumor_mut_file, 
    tumor_cna_file=tumor_cna_file, 
    tumor_exp_file=tumor_exp_file, 
    cell_line_mut_file=tumor_mut_file, 
    cell_line_cna_file=cell_line_cna_file, 
    cell_line_exp_file=cell_line_exp_file, 
    known_cancer_gene_weights_mut_file=known_cancer_gene_weights_mut_file, 
    known_cancer_gene_weights_cna_file=known_cancer_gene_weights_cna_file, 
    known_cancer_gene_weights_exp_file=known_cancer_gene_weights_exp_file, 
    cancer_specific_gene_weights_mut_file=cancer_specific_gene_weights_mut_file, 
    cancer_specific_gene_weights_cna_file=cancer_specific_gene_weights_cna_file, 
    cancer_specific_gene_weights_exp_file=cancer_specific_gene_weights_exp_file,
    distance_similarity_measures=c("generalized_jaccard", 
                                   "generalized_jaccard", 
                                   "weighted_correlation"))
  
  saved_output <- readRDS(system.file("test_output", "ov_comparison_result.rds", package="tumorcomparer"))
  
  expect_identical(comparison_result, saved_output)
})

test_that("run_shiny_app", {
  # TODO Add test
  expect_equal(TRUE, TRUE)
})


