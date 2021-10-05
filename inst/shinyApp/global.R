options(shiny.maxRequestSize = 24*1024^2) # Line exists in server, ui, global

# DEBUG ----
cat("MAX UPLOAD SIZE: ", getOption("shiny.maxRequestSize"), "\n")

# GENERAL ---- 
## cBioPortal 
cbioportal_mapping_file <- system.file(file.path("cbioportal_ccle", "ccle_cclp_cbioportal_mapping.txt"), package="tumorcomparer")
base_url_cbioportal <- '<a target="_blank" href="https://www.cbioportal.org/patient?studyId=ccle_broad_2019&caseId="'
cbioportal_mapping <- read.table(cbioportal_mapping_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cbioportal_mapping$CCLE_cBioPortal[!is.na(cbioportal_mapping$CCLE_cBioPortal)] <- 
  paste0(base_url_cbioportal, cbioportal_mapping$CCLE_cBioPortal[!is.na(cbioportal_mapping$CCLE_cBioPortal)], '">Link</a>')
colnames(cbioportal_mapping) <- c("Model_name", "TCGA_Type", "Info")

## Plotting
plot_title_prefix <-  "Mean Similarity to Tumors" 

plotlyModeBarButtonsToRemove <- c(
  "select2d", "sendDataToCloud", "pan2d", "resetScale2d", "hoverClosestCartesian", 
  "hoverCompareCartesian", "lasso2d", "zoomIn2d", "zoomOut2d", "toggleSpikelines"
)

# LOAD DATA ----
mtc_file <- system.file('extdata/mtc_results_20200331/mtc_results_20200331.rds', package="tumorcomparer")
#mtc_file <- system.file('extdata/mtc_results_20200331/mtc_results_20200331_no_factors.rds', package="tumorcomparer")
mtc_dataset <- readRDS(mtc_file)
precomputed_comparisons <- readRDS(system.file('extdata/precomputed_comparisons.rds', package="tumorcomparer"))
selected_geneset_comparions <- readRDS(system.file('extdata/selected_geneset_comparions.rds', package="tumorcomparer"))

mtc_dataset$Cell_Line_Name <- as.character(mtc_dataset$Cell_Line_Name)
#mtc_dataset$Cell_Line_Cancer_Type <- as.character(mtc_dataset$Cell_Line_Cancer_Type)
#mtc_dataset$Tumor_Cancer_Type <- as.character(mtc_dataset$Tumor_Cancer_Type)
mtc_dataset$Mean_Similarity_To_Tumors_AVG <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_AVG))[mtc_dataset$Mean_Similarity_To_Tumors_AVG]
mtc_dataset$Mean_Similarity_To_Tumors_MUT <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_MUT))[mtc_dataset$Mean_Similarity_To_Tumors_MUT]
mtc_dataset$Mean_Similarity_To_Tumors_CNA <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_CNA))[mtc_dataset$Mean_Similarity_To_Tumors_CNA] 
mtc_dataset$Mean_Similarity_To_Tumors_EXP <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_Tumors_EXP))[mtc_dataset$Mean_Similarity_To_Tumors_EXP]
mtc_dataset$AVGSIM_Zscores <- as.numeric(levels(mtc_dataset$AVGSIM_Zscores))[mtc_dataset$AVGSIM_Zscores]
mtc_dataset$MUTSIM_Zscores <- as.numeric(levels(mtc_dataset$MUTSIM_Zscores))[mtc_dataset$MUTSIM_Zscores]
mtc_dataset$CNASIM_Zscores <- as.numeric(levels(mtc_dataset$CNASIM_Zscores))[mtc_dataset$CNASIM_Zscores]
mtc_dataset$EXPSIM_Zscores <- as.numeric(levels(mtc_dataset$EXPSIM_Zscores))[mtc_dataset$EXPSIM_Zscores]
mtc_dataset$AVGSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$AVGSIM_Percentile_Ranks))[mtc_dataset$AVGSIM_Percentile_Ranks]
mtc_dataset$MUTSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$MUTSIM_Percentile_Ranks))[mtc_dataset$MUTSIM_Percentile_Ranks]
mtc_dataset$CNASIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$CNASIM_Percentile_Ranks))[mtc_dataset$CNASIM_Percentile_Ranks]
mtc_dataset$EXPSIM_Percentile_Ranks <- as.numeric(levels(mtc_dataset$EXPSIM_Percentile_Ranks))[mtc_dataset$EXPSIM_Percentile_Ranks]
mtc_dataset$Categorization <- as.character(mtc_dataset$Categorization)
mtc_dataset$AVGSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$AVGSIM_Zscores_wrt_Tumors))[mtc_dataset$AVGSIM_Zscores_wrt_Tumors]
mtc_dataset$MUTSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$MUTSIM_Zscores_wrt_Tumors))[mtc_dataset$MUTSIM_Zscores_wrt_Tumors]
mtc_dataset$CNASIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$CNASIM_Zscores_wrt_Tumors))[mtc_dataset$CNASIM_Zscores_wrt_Tumors]
mtc_dataset$EXPSIM_Zscores_wrt_Tumors <- as.numeric(levels(mtc_dataset$EXPSIM_Zscores_wrt_Tumors))[mtc_dataset$EXPSIM_Zscores_wrt_Tumors]
mtc_dataset$Mean_Similarity_To_All_Tumors_AVG <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_AVG))[mtc_dataset$Mean_Similarity_To_All_Tumors_AVG]
mtc_dataset$Mean_Similarity_To_All_Tumors_MUT <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_MUT))[mtc_dataset$Mean_Similarity_To_All_Tumors_MUT]
mtc_dataset$Mean_Similarity_To_All_Tumors_CNA <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_CNA))[mtc_dataset$Mean_Similarity_To_All_Tumors_CNA]
mtc_dataset$Mean_Similarity_To_All_Tumors_EXP <- as.numeric(levels(mtc_dataset$Mean_Similarity_To_All_Tumors_EXP))[mtc_dataset$Mean_Similarity_To_All_Tumors_EXP]
mtc_dataset$Average_Of_Percentile_Ranks <- as.numeric(levels(mtc_dataset$Average_Of_Percentile_Ranks))[mtc_dataset$Average_Of_Percentile_Ranks]
mtc_dataset$Rank_of_Average_Of_Percentile_Ranks <- as.numeric(levels(mtc_dataset$Rank_of_Average_Of_Percentile_Ranks))[mtc_dataset$Rank_of_Average_Of_Percentile_Ranks]

#saveRDS(mtc_dataset, 'inst/extdata/mtc_results_20200331/mtc_results_20200331_no_factors.rds')

# Column name mapping
mtc_selected_columns <- c(
  "Cell Line"="Cell_Line_Name", 
  #"Cell Line Type"="Cell_Line_Cancer_Type", 
  #"Tumor Type"="Tumor_Cancer_Type",
  "% Rank by Mutation"="MUTSIM_Percentile_Ranks", 
  "% Rank by Copy Number"="CNASIM_Percentile_Ranks", 
  "% Rank by Expression"="EXPSIM_Percentile_Ranks",
  "% Rank by Avg % Ranks"="Rank_of_Average_Of_Percentile_Ranks")

comparison_result_columns <- c(
  "Cell Line"="Cell_Line_Name", 
  #"Cell Line Type"="Cell_Line_Cancer_Type", 
  #"Tumor Type"="Tumor_Cancer_Type",
  "Mutation Score"="mut_score", 
  "Copy Number Score"="cna_score", 
  "Expression Score"="exp_score",
  "Combined Score"="combined_score") 

tcgaTypes <- c(
  "Adrenocortical Carcinoma"="ACC",
  "Bladder Urothelial Carcinoma"="BLCA",
  "Breast Invasive Carcinoma"="BRCA",
  "Cervical Squamous Cell Carcinoma"="CESC",
  "Cholangiocarcinoma"="CHOL",
  "Colon Adenocarcinoma"="COAD",
  "Diffuse Large B-Cell Lymphoma"="DLBC",
  "Esophageal Adenocarcinoma"="ESCA",
  "Glioblastoma Multiforme"="GBM", 
  "Head and Neck Squamous Cell Carcinoma"="HNSC",
  "Kidney Renal Clear Cell Carcinoma"="KIRC",
  "Acute Myeloid Leukemia"="LAML",
  "Low-Grade Glioma"="LGG",
  "Liver Hepatocellular Carcinoma"="LIHC",
  "Lung adenocarcinoma"="LUAD",   
  "Lung squamous cell"="LUSC", 
  "Mesothelioma"="MESO",
  "Ovarian"="OV",
  "Pancreatic Adenocarcinoma"="PAAD",
  "Prostate Adenocarcinoma"="PRAD",
  "Rectal Adenocarcinoma"="READ",
  "Cutaneous Melanoma"="SKCM",
  "Stomach Adenocarcinoma"="STAD",
  "Thyroid Cancer"="THCA",
  "Endometrial Carcinoma"="UCEC"
)
tcgaTypes <- tcgaTypes[tcgaTypes %in% as.character(unique(mtc_dataset$Tumor_Cancer_Type))]


genesets <- c("Most Variable Genes", "Cell Cycle", "HIPPO", "MYC", "NOTCH", "PI3K", "RTK RAS", "TGF-Beta", "WNT")