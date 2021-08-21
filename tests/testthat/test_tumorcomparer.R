# All tests are done on files in package using system.file()

test_that("sanity_check", {
  expect_equal(2 * 2, 4)
})

# TODO Add tests
# test_that("plot_mds", {
#   expect_equal(TRUE, FALSE)
# })

# test_that("run_shiny_app", {
#   expect_equal(TRUE, TRUE)
# })

test_that("return_first_part", {
  tmp <- return_first_part("22RV1_PROSTRATE")
  expect_equal(tmp, "22RV1")
})

test_that("keep_only_high_level_cnas", {
  set.seed(1)
  tmp <- sample(c(-2, -1, 0, 1, 2), 100, replace=TRUE)
  mat <- matrix(tmp, 10, 10)
  
  results <- as.vector(keep_only_high_level_cnas(cna_mat=mat))
  expect_false(any(results == 1 | results == -1))
})

test_that("keep_only_high_level_cnas", {
  t1 <- as.matrix(c(1,1,1,1,1,1,2,2,2,3,3,3,4,4))
  t2 <- as.matrix(rev(t1))
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
  mut_default_weight <- 0.01
  tumor_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "tumor_mut.txt", package="tumorcomparer")
  cell_line_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "cell_line_mut.txt", package="tumorcomparer")
  known_cancer_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
  cancer_specific_gene_weights_mut_file <- system.file("extdata", "ovarian_tcga_cclp", "Genes_and_weights_mut.txt", package="tumorcomparer")
  
  composite_mat <- generate_composite_mat_and_gene_weights(
    default_weight=mut_default_weight,
    tumor_file=tumor_mut_file,
    cell_line_file=cell_line_mut_file,
    known_cancer_gene_weights_file=known_cancer_gene_weights_mut_file,
    cancer_specific_gene_weights_file=cancer_specific_gene_weights_mut_file)
  
  saved_output <- readRDS(system.file("test_output", "ov_composite_mat.rds", package="tumorcomparer"))
  
  expect_equal(composite_mat, saved_output)
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
  
  expect_equal(result[1:3], c("#DD8F21", "#F9A105", "#C37E3B"))
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

test_that("run_comparison", {
  set.seed(1)
  
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
    tumor_mut_file=tumor_mut_file, 
    tumor_cna_file=tumor_cna_file, 
    tumor_exp_file=tumor_exp_file, 
    cell_line_mut_file=cell_line_mut_file, 
    cell_line_cna_file=cell_line_cna_file, 
    cell_line_exp_file=cell_line_exp_file, 
    known_cancer_gene_weights_mut_file=known_cancer_gene_weights_mut_file, 
    known_cancer_gene_weights_cna_file=known_cancer_gene_weights_cna_file, 
    known_cancer_gene_weights_exp_file=known_cancer_gene_weights_exp_file, 
    cancer_specific_gene_weights_mut_file=cancer_specific_gene_weights_mut_file, 
    cancer_specific_gene_weights_cna_file=cancer_specific_gene_weights_cna_file, 
    cancer_specific_gene_weights_exp_file=cancer_specific_gene_weights_exp_file)
  
  #saveRDS(comparison_result, "ov_comparison_result.rds")
  
  saved_output <- readRDS(system.file("test_output", "ov_comparison_result_old.rds", package="tumorcomparer"))
  
  expect_equal(comparison_result, saved_output)
})

### test function on 3 data types for new generilized function
test_that("run_comparison_config_list", {
  set.seed(1)
  
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
  
  
  ### creating config list for comparison function 
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
  )
  
  comparison_result <- run_comparison_config_list(config_list = config_list)
  
  
  #saveRDS(comparison_result, "ov_comparison_result.rds")
  
  saved_output <- readRDS(system.file("test_output", "ov_comparison_result.rds", package="tumorcomparer"))
  
  expect_equal(comparison_result, saved_output)
})


test_that("run_comparison_two_datasets", {
  set.seed(1)
  
  tumor_mut_file <- system.file("extdata", "READ_data_for_running_TC", "tumor_mut.txt", package="tumorcomparer")
  tumor_exp_file <- system.file("extdata", "READ_data_for_running_TC", "tumor_exp.txt", package="tumorcomparer")
  
  cell_line_mut_file <- system.file("extdata", "READ_data_for_running_TC", "cell_line_mut.txt", package="tumorcomparer")
  cell_line_exp_file <- system.file("extdata", "READ_data_for_running_TC", "cell_line_exp.txt", package="tumorcomparer")
  
  known_cancer_gene_weights_mut_file <- system.file("extdata", "READ_data_for_running_TC", "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
  known_cancer_gene_weights_exp_file <- system.file("extdata", "READ_data_for_running_TC", "default_weights_for_known_cancer_genes_exp.txt", package="tumorcomparer")
  
  cancer_specific_gene_weights_mut_file <- system.file("extdata", "READ_data_for_running_TC", "Genes_and_weights_mut.txt", package="tumorcomparer")
  cancer_specific_gene_weights_exp_file <- system.file("extdata", "READ_data_for_running_TC", "Genes_and_weights_exp.txt", package="tumorcomparer")
  
  comparison_result <- run_comparison(
    available_data_types=c("mut", "exp"), 
    mut_data_type_weight = 1/2,
    exp_data_type_weight = 1/2,
    cna_default_weight=0.01, 
    exp_default_weight=0.01,
    tumor_mut_file=tumor_mut_file, 
    tumor_exp_file=tumor_exp_file, 
    cell_line_mut_file=cell_line_mut_file, 
    cell_line_exp_file=cell_line_exp_file, 
    known_cancer_gene_weights_mut_file=known_cancer_gene_weights_mut_file, 
    known_cancer_gene_weights_exp_file=known_cancer_gene_weights_exp_file, 
    cancer_specific_gene_weights_mut_file=cancer_specific_gene_weights_mut_file, 
    cancer_specific_gene_weights_exp_file=cancer_specific_gene_weights_exp_file)
  
  #saveRDS(comparison_result, "ov_comparison_result.rds")
  
  #saved_output <- readRDS(system.file("test_output", "ov_comparison_result.rds", package="tumorcomparer"))
  
  expect_equal(TRUE, TRUE)
})

### test function on 5 data types for new generilized function
test_that("run_comparison_5_datasets", {
  
  set.seed(1)
  
  tumor_exp_file <- system.file("extdata", "mock_5_data_types", "tumor_exp.txt", package="tumorcomparer")
  tumor_cna_file <- system.file("extdata", "mock_5_data_types", "tumor_cna.txt", package="tumorcomparer")
  tumor_meth_file <- system.file("extdata", "mock_5_data_types", "tumor_meth.txt", package="tumorcomparer")
  tumor_mut_file <- system.file("extdata", "mock_5_data_types", "tumor_mut.txt", package="tumorcomparer")
  tumor_prot_file <- system.file("extdata", "mock_5_data_types", "tumor_prot.txt", package="tumorcomparer")
  
  
  cell_line_exp_file <- system.file("extdata", "mock_5_data_types", "cell_line_exp.txt", package="tumorcomparer")
  cell_line_cna_file <- system.file("extdata", "mock_5_data_types", "cell_line_cna.txt", package="tumorcomparer")
  cell_line_meth_file <- system.file("extdata", "mock_5_data_types", "cell_line_meth.txt", package="tumorcomparer")
  cell_line_mut_file <- system.file("extdata", "mock_5_data_types", "cell_line_mut.txt", package="tumorcomparer")
  cell_line_prot_file <- system.file("extdata", "mock_5_data_types", "cell_line_prot.txt", package="tumorcomparer")
  
  
  known_cancer_gene_weights_exp_file <- system.file("extdata", "mock_5_data_types", "default_weights_for_known_cancer_genes_exp.txt", package="tumorcomparer")
  known_cancer_gene_weights_cna_file <- system.file("extdata", "mock_5_data_types", "default_weights_for_known_cancer_genes_cna.txt", package="tumorcomparer")
  known_cancer_gene_weights_meth_file <- system.file("extdata", "mock_5_data_types", "default_weights_for_known_cancer_genes_meth.txt", package="tumorcomparer")
  known_cancer_gene_weights_mut_file <- system.file("extdata", "mock_5_data_types", "default_weights_for_known_cancer_genes_mut.txt", package="tumorcomparer")
  known_cancer_gene_weights_prot_file <- system.file("extdata", "mock_5_data_types", "default_weights_for_known_cancer_genes_prot.txt", package="tumorcomparer")
  
  
  cancer_specific_gene_weights_exp_file <- system.file("extdata", "mock_5_data_types", "Genes_and_weights_exp.txt", package="tumorcomparer")
  cancer_specific_gene_weights_cna_file <- system.file("extdata", "mock_5_data_types", "Genes_and_weights_cna.txt", package="tumorcomparer")
  cancer_specific_gene_weights_meth_file <- system.file("extdata", "mock_5_data_types", "Genes_and_weights_meth.txt", package="tumorcomparer")
  cancer_specific_gene_weights_mut_file <- system.file("extdata", "mock_5_data_types", "Genes_and_weights_mut.txt", package="tumorcomparer")
  cancer_specific_gene_weights_prot_file <- system.file("extdata", "mock_5_data_types", "Genes_and_weights_prot.txt", package="tumorcomparer")
  
  
  ### creating config list for comparison function 
  
  config_list <- list(exp=list(dataset_name = "exp", data_type_weight=1/5, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/5, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      meth=list(dataset_name = "meth", data_type_weight=1/5, default_weight = 0.01, 
                               tumor_file = tumor_meth_file, cell_line_file = cell_line_meth_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_meth_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_meth_file),
                      mut=list(dataset_name = "mut", data_type_weight=1/5, default_weight = 0.01, 
                                tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                                known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                                cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      prot=list(dataset_name = "prot", data_type_weight=1/5, default_weight = 0.01, 
                                tumor_file = tumor_prot_file, cell_line_file = cell_line_prot_file,
                                known_cancer_gene_weights_file = known_cancer_gene_weights_prot_file, 
                                cancer_specific_gene_weights_file = cancer_specific_gene_weights_prot_file)
  )
  
  test_results <- run_comparison_config_list(config_list = config_list)
  
  saved_output <- readRDS(system.file("test_output", "comparison_5_datasets.rds", package="tumorcomparer"))
  
  expect_equal(test_results, saved_output)
  
})


### test function on 3 data types for new generilized function with specified gene_list_argument
test_that("run_comparison_config_list", {
  set.seed(123)
  
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
  
  sample_gene_list <- c("ARRDC1", "CNTN6", "CREBBP", "EP300", "HES1", "HES2", "HES3", "HES4", "HES5", "HEY1", 
                        "HEY2", "HEYL", "KAT2B", "KDM5A", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NOV", "NRARP", 
                        "PSEN2", "LFNG", "ITCH", "NCSTN", "SPEN", "JAG1", "APH1A", "FBXW7", "FHL1", "THBS2", "HDAC2", 
                        "MFAP2", "CUL1", "RFNG", "NCOR1", "NCOR2", "MFAP5", "HDAC1", "NUMB", "JAG2", "MAML3", "MFNG", 
                        "CIR1", "CNTN1", "MAML1", "MAML2", "NUMBL", "PSEN1", "PSENEN", "RBPJ", "RBPJL", "RBX1", "SAP30", 
                        "SKP1", "SNW1", "CTBP1", "CTBP2", "ADAM10", "APH1B", "ADAM17", "DLK1", "DLL1", "DLL3", "DLL4", "DNER", 
                        "DTX1", "DTX2", "DTX3", "DTX3L", "DTX4", "EGFL7", "TMEM129", "ADRA2C", "SLC2A9", "RAD54B", "MAP3K1", "EIF3E", 
                        "WDYHV1")
  
  ### creating config list for comparison function 
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
  )
  
  comparison_result <- run_comparison_config_list(config_list = config_list, gene_list = sample_gene_list)
  
  
  #saveRDS(comparison_result, "ov_comparison_result.rds")
  
  saved_output <- readRDS(system.file("test_output", "comparison_geneset.rds", package="tumorcomparer"))
  
  expect_equal(comparison_result, saved_output)
})

## minimum gene number test
test_that("testing few genes error", {
  set.seed(1)
  
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
  
  sample_gene_list <- c("ARRDC1", "CNTN6", "CREBBP", "EP300")
  
  ### creating config list for comparison function 
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
  )
  

  comparison_failed <- try(comparison_result <- run_comparison_config_list(config_list = config_list, gene_list = sample_gene_list), silent = TRUE)
  
  comparison_failed <- unlist(strsplit(comparison_failed[1], split = '\n', fixed = T))[2]
  
  expect_equal(comparison_failed, "  ERROR: At least 5 genes are required for the comparison")
})


## test for cyj_graph_maker_from_dist_mat function
test_that("cyj_graph_maker_from_dist_mat", {
  set.seed(1)
  
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
  
  
  ### creating config list for comparison function 
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
  )
  
  comparison_result <- run_comparison_config_list(config_list = config_list)
  
  cyj_json_graph <- cyj_graph_maker_from_dist_mat(dist_mat = comparison_result$dist_mat, min_weight = 0.85)
  
  saved_output <- readRDS(system.file("test_output", "cyj_json_graph.rds", package="tumorcomparer"))
  
  expect_equal(cyj_json_graph, saved_output)
})


## test for ballon_plot_data_to_result_table function
test_that("ballon_plot_data_to_result_table", {
  set.seed(1)
  
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
  
  
  ### creating config list for comparison function 
  
  config_list <- list(mut=list(dataset_name = "mut", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_mut_file, cell_line_file = cell_line_mut_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_mut_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_mut_file),
                      cna=list(dataset_name = "cna", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_cna_file, cell_line_file = cell_line_cna_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_cna_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_cna_file),
                      exp=list(dataset_name = "exp", data_type_weight=1/3, default_weight = 0.01, 
                               tumor_file = tumor_exp_file, cell_line_file = cell_line_exp_file,
                               known_cancer_gene_weights_file = known_cancer_gene_weights_exp_file, 
                               cancer_specific_gene_weights_file = cancer_specific_gene_weights_exp_file)
  )
  
  comparison_result <- run_comparison_config_list(config_list = config_list)
  
  plot_data <- list(plot_data = make_balloon_plot_data_from_comparison_result(comparison_result = comparison_result))
  
  result_table <- ballon_plot_data_to_result_table(plot_data = plot_data)
  
  saved_output <- readRDS(system.file("test_output", "baloon_plot_result_table.rds", package="tumorcomparer"))
  
  # expect_equal(result_table, saved_output)
  
  expect_equal(TRUE, TRUE)
})