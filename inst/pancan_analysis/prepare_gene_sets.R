# Set working directory
setwd("~/Google Drive/tumorcomparer/pancan_data")

# Prepares gene sets 
gene_sets <- list()

# Vogelstein lists 
vogelstein <- read.table("Vogelstein_Drivers_Mutn_and_CNAs.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)

tsg_idx <- which(vogelstein$Category == "TSG")
og_idx <- which(vogelstein$Category == "Oncogene")

gene_sets[["Tumor Suppressor"]] <- vogelstein[tsg_idx, "Gene"]
gene_sets[["Oncogene"]] <- vogelstein[og_idx, "Gene"]

# Kinases 
kinases <- read.table("Human_Kinase_List.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
gene_sets[["Kinase"]] <- kinases[, "Official.Symbol"]

# Phosphatases
phosphatases <- read.table("PhosphoClass_2017-05-19_21-30.csv", sep=";", header=TRUE, stringsAsFactors = FALSE)
gene_sets[["Phosphatase"]] <- trimws(phosphatases[, "Protein.Name"])

# Convert to data.frame
df <- data.frame(gene=character(0), category=character(0))

for(i in 1:length(gene_sets)) {
  t1 <- gene_sets[[names(gene_sets[i])]]
  t2 <- data.frame(gene=t1, category=names(gene_sets[i]))
  df <- rbind(df, t2)
}

# Export
datestamp <- format(Sys.time(), "%m-%d-%Y")
saveRDS(df, paste0("gene_sets/gene_sets_", datestamp, ".rds"))
