library(tidyverse)

# SET WORKING DIRECTORY ----
setwd("~/Google\ Drive/tumorcomparer/pancan_data/")
stopifnot(file.exists("COSMIC_v81/CosmicMutantExport.tsv"))

# Parameters 
repeat_cutoff <- 10 

# Read data
cosmic <- read_tsv("COSMIC_v81/CosmicMutantExport.tsv", 
                   col_names = TRUE, na = "")

# Extract needed data
tmp <- cosmic[, c("Gene name", "Mutation AA")]

# Add column for count of entries is original file 
t2 <- summarise(group_by(tmp, `Gene name`, `Mutation AA`), length(`Mutation AA`))
colnames(t2)[3] <- "COSMIC_StudyCount"

# Validation check
t3 <- t2[which(t2$`Gene name` == "BRAF" & t2$`Mutation AA` == "p.V600E"), ]
stopifnot(t3$COSMIC_StudyCount > 45000)

filtered_cosmic_mutations <- t2[t2$COSMIC_StudyCount >= repeat_cutoff,]
no_question_mark_idx <- which(!grepl("\\?", filtered_cosmic_mutations$`Mutation AA`))

filtered_cosmic_mutations <- filtered_cosmic_mutations[no_question_mark_idx, ]
nrow(filtered_cosmic_mutations)

saveRDS(filtered_cosmic_mutations, "filtered_cosmic_mutations.rds")
