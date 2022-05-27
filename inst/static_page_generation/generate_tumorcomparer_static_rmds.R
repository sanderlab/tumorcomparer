library(brew)
library(rmarkdown)
library(stringr)
library(xml2)
library(magrittr)

work_dir <- "inst/static_page_generation/"
setwd(work_dir)

source("https://raw.githubusercontent.com/cannin/cbioportal_to_complexheatmap_oncoprint/master/transform_cbioportal_to_complexheatmap_oncoprint.R")

# Parameters 
labels <- c("mRNA High", "mRNA Low", "Protein High",
            "Protein Low", "Missense Mutation (putative passenger)", "Missense Mutation (putative driver)", 
            "Inframe Mutation (putative passenger)", "Inframe Mutation (putative driver)", "Truncating mutation (putative passenger)",
            "Truncating mutation (putative driver)", "Amplification", "Deep Deletion", 
            "homdel_rec", "amp_rec", "sv", 
            "sv_rec", "splice")
names(labels) <- c("mrna_hi", "mrna_lo", "prot_hi", 
                   "prot_lo", "mut_mis_put_pass", "mut_mis_put_driv", 
                   "mut_inframe_put_pass", "mut_inframe_put_driv", "mut_trun_put_pass", 
                   "mut_trun_put_driv", "amp", "deep_del", 
                   "homdel_rec", "amp_rec", "sv", 
                   "sv_rec", "splice")

oncoprint_caption <- paste(sapply(1:length(labels), function(i) {
  paste0(names(labels[i]), " = ", labels[i])
}), collapse = "; ")

# Function to generate report and render output
create.report <- function(reports, prepend=NULL) {
  
  # Set parameters 
  file_prefix <- reports$tcga_abbrev

  tcga_title <- str_to_title(reports$tcga_title)
  tcga_abbrev <- reports$tcga_abbrev
  
  tmp_oncoprint_caption <- paste0(c("Comparing", tolower(tcga_title), "cell lines to tumor samples."), collapse=" ")
  chunk_oncoprint_caption <- paste0(tmp_oncoprint_caption, " Legend Abbreviations: ", oncoprint_caption)
  chunk_oncoprint_name <- tolower(gsub(" ", "_", tmp_oncoprint_caption))
  
  rmd_file <- paste0(file_prefix, "_cclp_tcga_tumorcomparer.Rmd")
  
  #DEBUG
  cat(".Rmd: ", rmd_file, "\n")
  
  brew("tumorcomparer_template.brew", rmd_file)
  
  # Generate report
  #render(rmd_file)		
}

files <- dir(".", "^[a-z]*_cclp_tcga_tumorcomparer.Rmd")
file.remove(files)

ccle_tcga_mapping <- read.table("ccle_tcga_mapping.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

reports <- data.frame(tcga_abbrev=tolower(ccle_tcga_mapping$tcga), tcga_title=ccle_tcga_mapping$name, stringsAsFactors=FALSE)

# GENERATE REPORTS ----
for(i in 1:nrow(reports)) {
  create.report(reports[i,])
}

# RENDER SITE ----
rmarkdown::render_site()

# CREATE SITEMAP ----
output_dir <- "_static"
base_url <- "http://projects.sanderlab.org/tumorcomparer/static/"
file_list <- list.files(output_dir, "html")
#writeLines(paste0(base_url, file_list), file.path(output_dir, "sitemap.txt"))

# Sitemap XML Format Description: https://www.sitemaps.org/protocol.html
doc <- xml_new_root("urlset", xmlns="http://www.sitemaps.org/schemas/sitemap/0.9") 

for(page in file_list) {
  url <- xml_new_document() %>% xml_add_child("url")
  
  loc <- xml_new_document() %>% xml_add_child("loc")
  xml_text(loc) <- paste0(base_url, page)
  
  lastmod <- xml_new_document() %>% xml_add_child("lastmod")
  xml_text(lastmod) <- Sys.Date() %>% as.character
  
  changefreq <- xml_new_document() %>% xml_add_child("changefreq")
  xml_text(changefreq) <- "monthly"
  
  xml_add_child(url, loc)
  xml_add_child(url, lastmod)
  xml_add_child(url, changefreq)
  
  xml_add_child(doc, url)  
}

invisible(write_xml(doc, file.path(output_dir, "sitemap.xml")))



#source("generate_tumorcomparer_static_rmds.R")

# CLEANUP ----
file.remove(list.files(pattern="*_tcga_tumorcomparer.[Rmd|html]", full.names=TRUE))
#file.remove(list.files(pattern="*_tcga_tumorcomparer", full.names=TRUE))

