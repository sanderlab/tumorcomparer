library(magrittr)

tmp <- read.table("inst/cbioportal_ccle/cellline_ccle_broad_clinical_data.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
tmp$tmp_ccle <- tmp$Patient.ID %>% tolower %>% gsub("[^[:alnum:]]", "", .)

dat <- read.table("inst/cbioportal_ccle/CCLP_ID_and_Cancer_Type.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dat$tmp_cclp <- dat$Model_name %>% tolower %>% gsub("[^[:alnum:]]", "", .)

x1 <- merge(dat, tmp, by.x="tmp_cclp", by.y="tmp_ccle", all.x=TRUE)
x1 <- x1[, c("Model_name", "TCGA_Type", "Sample.ID")]
colnames(x1) <- c("Model_name", "TCGA_Type", "CCLE_cBioPortal")

write.table(x1, file="inst/cbioportal_ccle/ccle_cclp_cbioportal_mapping.txt", sep="\t", row.names=FALSE, quote=FALSE)
