#--------------------------------------------------------------------------------------------------
# LOAD AND CHECK DATA
#--------------------------------------------------------------------------------------------------
library(tidyverse)

workDir <- "~/default/workspaceNotSynced/annovar/tumorcomparer_annovar_inputs/gdsc/"
gdscFile <- file.path(workDir, "CosmicCLP_MutantExport.tsv")
fullMut <- read.table(file=gdscFile, header=TRUE, sep="\t", comment.char="", quote="", stringsAsFactors=FALSE)

fullMut <- filter(fullMut, fullMut$Mutation.Description != "Substitution - coding silent")

# Additional column with actual gene names (original Gene.name column sometimes has
# GENE_TRANSCRIPT_ID names).
fullMut <- rename(fullMut, GENE_TRANSCRIPT = Gene.name, CHROM_POS = Mutation.genome.position)
fullMut$GENE <- vapply(fullMut$GENE_TRANSCRIPT, function(x){
  stringr::str_split(x, "_")[[1]][1]
}, character(1))
fullMut <- fullMut[, c(39, 1:38)]

#--------------------------------------------------------------------------------------------------
# CONSTRUCT MISSENSE MUTATION DATA TABLE, ANNOVAR INPUT
#--------------------------------------------------------------------------------------------------
#-----[helper functions]------------------------------------------------------
getPairedBase <- function(x){
  pairedBase <- NA
  if (x == "A"){
    pairedBase <- "T"
  } else if (x == "T"){
    pairedBase <- "A"
  } else if (x == "G"){
    pairedBase <- "C"
  } else if (x == "C"){
    pairedBase <- "G"
  }
  return(pairedBase)
}

# Assumes MISSENSE mutations with format "c.3347G>T"
getRefObs <- function(varStr, strand){
  results <- list()
  results$ref <- NA
  results$obs <- NA

  tmp <- stringr::str_split(varStr, pattern = "\\d+")[[1]][2]
  if (nchar(tmp) == 3){
    tmp <- stringr::str_extract(tmp, pattern = "[ATGC]>[ATGC]$")
  } else{
    tmp <- NA
  }

  if (!is.na(tmp)){
    tmp <- stringr::str_split(tmp, pattern = ">")[[1]]
    if (strand == "+"){
      results$ref <- tmp[1]
      results$obs <- tmp[2]
    } else if (strand == "-"){
      results$ref <- getPairedBase(tmp[1])
      results$obs <- getPairedBase(tmp[2])
    }
  }

  return(results)
}

# Assumes single base (MISSENSE) mutation chromosomal locations
# with format "7:140753336-140753336"
getChromLoc <- function(chrStr){
  results <- list()
  results$chr <- NA
  results$pos <- NA

  tmp <- stringr::str_split(chrStr, pattern = ":")[[1]]
  chr <- tmp[1]
  posStr <- tmp[2]
  pos <- unique(stringr::str_split(posStr, pattern = "-")[[1]][1])
  if (length(pos) == 1){
    results$chr <- chr
    results$pos <- as.integer(pos)
  }

  return(results)
}

#------------------------------------------------------------------------------

missMut <- filter(fullMut, Mutation.Description == "Substitution - Missense")
uniqMissMut <- select(missMut, GENE, Mutation.CDS, CHROM_POS, strand, GRCh) %>% unique()
uniqMissMut$CHROM <- NA
uniqMissMut$START <- NA
uniqMissMut$END <- NA
uniqMissMut$REF <- NA
uniqMissMut$OBS <- NA
uniqMissMut$ID <- NA

N <- nrow(uniqMissMut)
for (i in (1:N)){
  varStr <- uniqMissMut[i, "Mutation.CDS"]
  chrStr <- uniqMissMut[i, "CHROM_POS"]
  strand <- uniqMissMut[i, "strand"]

  chrLoc <- getChromLoc(chrStr)
  if (!any(is.na(chrLoc))){
    uniqMissMut[i, "CHROM"] <- chrLoc$chr
    uniqMissMut[i, "START"] <- chrLoc$pos
    uniqMissMut[i, "END"]   <- chrLoc$pos
  }

  refObs <- getRefObs(varStr, strand)
  if (!any(is.na(refObs))){
    uniqMissMut[i, "REF"] <- refObs$ref
    uniqMissMut[i, "OBS"] <- refObs$obs
  }

  uniqMissMut[i, "ID"] <- stringr::str_c(
    uniqMissMut[i, c("GENE", "CHROM", "START", "END", "REF", "OBS")], collapse = "_")
  if ((i %% 100) == 0){
    cat(paste0(signif((i/N)*100, 3), "% completed."), sep = "\n")
  }
}

save(uniqMissMut, file = file.path(workDir, "uniqMissMut.RData"))

#--------------------------------------------------------------------------------------------------
# WRITE ANNOVAR INPUT FILE (MISSENSE MUTATIONS ONLY)
#--------------------------------------------------------------------------------------------------
load(file.path(workDir, "uniqMissMut.RData"))

tmp <- filter(uniqMissMut, !is.na(GRCh))
stopifnot(all(tmp$GRCh == 37))
tmp <- select(tmp, CHROM, START, END, REF, OBS, ID)

# Make sure all columns are free of NA values.
stopifnot(all(unlist(vapply(tmp, function(x) {!any(is.na(x))}, logical(1)))))

write_tsv(tmp, path = file.path(workDir, "gdsc_msmut_annovar_input.txt"), col_names=FALSE)

# ----[ANNOVAR COMMAND (concatenate to single line)]-------------------------------------
# perl table_annovar.pl 
# ~/REPOS/gdscdatadec15/inst/extdata/cosmic_v79/gdsc_missense_mut_annovar_input.txt  
# humandb/ -buildver hg38 -out gdsc_msmut_annovar -remove 
# -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,ljb26_all 
# -operation g,r,r,f,f,f -nastring -