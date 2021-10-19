library(tumorcomparer)
library(data.table)
library(simpleRCache)
library(magrittr)

#' Starts a stopwatch timer to measure performance
#' 
#' @param gcFirst a boolean whether to run garbage collection before starting watch
#' @param type the time of time to keep track of 
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

#' Stops a stopwatch timer to measure performance 
#' 
toc <- function() {
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  diff <- toc - tic
  print(diff)
  invisible(toc)
  return(diff)
}

data("tcga_pancan_pathway_genes")

# PARAMETERS ----
## Caching
tc_geneset_comparison_cached <- addMemoization(tc_geneset_comparison)

## SETUP CACHE
Sys.setenv("DEBUG_SIMPLERCACHE"="TRUE")
cache_dir <- ifelse(
  Sys.getenv("CACHE_DIR") != "", 
  Sys.getenv("CACHE_DIR"), 
  'cache_tc')
setCacheRootPath(cache_dir)

tc_dataset_dir <- "~/default/workspaceNotSynced/tumorcomparer_data/TC_Data_PanCancer_March2021/"

most_variable_genes_precomputed_path <- system.file("extdata", 
                                                    "mtc_results_20200331", 
                                                    "mtc_results_20200331_no_factors.rds", 
                                                    package = "tumorcomparer")
most_variable_genes_precomputed_results <- readRDS(most_variable_genes_precomputed_path)
avail_cancer_types <- as.character(unique(most_variable_genes_precomputed_results$Tumor_Cancer_Type))

gene_lists <- names(tcga_pancan_pathway_genes)
drop_entries <- c("NRF2")
gene_lists <- gene_lists[!(gene_lists %in% drop_entries)]

remove_tmp_files <- TRUE
remove_errored_dataset_comparisons <- TRUE

entries <- expand.grid(gene_lists, avail_cancer_types, stringsAsFactors=FALSE)
entries$time <- NA
entries$file <- NA
colnames(entries) <- c("gene_list", "cancer_type", "time", "file")
#entries <- entries[entries$gene_list == "TP53", ]

for(i in 1:nrow(entries)) {
#for(i in 1:3) {
  cancer_type <- entries$cancer_type[i]
  gene_list <- entries$gene_list[i]
  
  # gene_list <- "TP53"; cancer_type <- "BLCA"

  #cancer_type <- "LIHC"
  genes <- tcga_pancan_pathway_genes[[gene_list]]
  
  tryCatch({
    tic()
    
    file_prefix <- tolower(paste0(gsub("[^[:alnum:]]", "_", gene_list), "_", cancer_type))
    Sys.setenv("PREFIX_SIMPLERCACHE"=file_prefix)
    
    results_hash <- capture.output({ 
      comparison_result <- tc_geneset_comparison_cached(
        gene_list=genes, 
        cancer_type=cancer_type, 
        tc_dataset_dir=tc_dataset_dir, 
        remove_errored_dataset_comparisons=remove_errored_dataset_comparisons, 
        remove_tmp_files=remove_tmp_files,
        verbose=TRUE)      
      
    })
    
    results_file <- results_hash[2] %>% trimws %>% sub("DEBUG: Cache used:  cache_tc/", "", .)
    
    diff_time <- toc()
    
    entries$time[i] <- diff_time
    entries$file[i] <- results_file
    status_str <- paste0("I:", i, " G: ", gene_list, "; C: ", cancer_type, "; T: ", round(diff_time, 2), " F: ", results_file, "\n")
    cat(status_str)
  }, error=function(e) {
    status_str <- paste0("ERROR: I:", i, " G: ", gene_list, " C: ", cancer_type, "; T: NA\n")
    cat(status_str)
    cat("ERROR: ", e$message, "\n")
  })
}

stopifnot(!any(is.na(entries$time)))

# Process Pre-computed Comparisons
tmp_cancer_types <- unique(entries$cancer_type)
tmp_gene_lists <- unique(entries$gene_list)

precomputed_comparisons <- lapply(tmp_cancer_types, function(x) {
  lapply(tmp_gene_lists, function(y) {

    print(x)
    print(y)
    
    cur_file <- entries$file[entries$cancer_type == x & entries$gene_list == y]
    
    comparison_result <- readRDS(file.path(cache_dir, cur_file))

    plot_data <- tumorcomparer::make_balloon_plot_data_from_comparison_result(comparison_result)

    list(compared_data_types = comparison_result$calculated_data_types, plot_data = plot_data)
  })
})

names(precomputed_comparisons) <- tmp_cancer_types

precomputed_comparisons <- lapply(precomputed_comparisons, function(x) {
  names(x) <- tmp_gene_lists
  x
})

saveRDS(precomputed_comparisons, file="inst/extdata/precomputed_geneset_comparisons/precomputed_comparisons_20211019.rds")

# Make selected_geneset_comparisons
tmp_gene_lists <- c("Most Variable Genes", tmp_gene_lists)

selected_geneset_comparisons <- lapply(tmp_cancer_types, function(x) {
  tmp_gene_lists
})

names(selected_geneset_comparisons) <- tmp_cancer_types

saveRDS(selected_geneset_comparisons, file="inst/extdata/precomputed_geneset_comparisons/selected_geneset_comparisons_20211019.rds")



