# GENERAL 

plot_title_prefix <-  "Mean Similarity to Tumors" 

plotlyModeBarButtonsToRemove <- c(
  "select2d", "sendDataToCloud", "pan2d", "resetScale2d", "hoverClosestCartesian", 
  "hoverCompareCartesian", "lasso2d", "zoomIn2d", "zoomOut2d", "toggleSpikelines"
)

# LOAD DATA ----
mtc_file <- system.file('extdata/mtc_results_20200331/mtc_results_20200331.rds', package="tumorcomparer")
mtc_dataset <- readRDS(mtc_file)
mtc_selected_columns <- c(
  "Cell Line"="Cell_Line_Name", 
  #"Cell Line Type"="Cell_Line_Cancer_Type", 
  #"Tumor Type"="Tumor_Cancer_Type",
  "% Rank by Mutation"="MUTSIM_Percentile_Ranks", 
  "% Rank by Copy Number"="CNASIM_Percentile_Ranks", 
  "% Rank by Expression"="EXPSIM_Percentile_Ranks",
  "% Rank by Avg % Ranks"="Rank_of_Average_Of_Percentile_Ranks")

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

