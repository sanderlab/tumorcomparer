# LOAD DATA
#lapply(dir(file.path("www", "db"), pattern="RData", recursive=TRUE, full.names=TRUE), load, .GlobalEnv)

#tcgaTypes <- c("Ovarian (OV)"="OV", "Lung squamous cell (LUSC)"="LUSC", "Lung adenocarcinoma (LUAD)"="LUAD", "Glioblastoma multiforme (GBM)"="GBM", "Colorectal adenocarcinoma (COADREAD)"="COADREAD", "Breast invasive carcinoma (BRCA)"="BRCA")
tcgaTypes <- c("Lung squamous cell"="LUSC", "Lung adenocarcinoma"="LUAD", "Ovarian"="OV", "Glioblastoma multiforme"="GBM", "Colorectal adenocarcinoma"="COADREAD", "Breast invasive carcinoma"="BRCA")

categorizations <- c("Great", "Good", "Moderately Good", "Poor", "Outlier")