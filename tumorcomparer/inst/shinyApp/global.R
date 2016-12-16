# LOAD DATA
lapply(dir(file.path("www", "db"), pattern="RData", recursive=TRUE, full.names=TRUE), load, .GlobalEnv)

# DEBUGGING
# lapply(dir(
#     file.path("inst", "shinyApp", "www", "db"),
#     pattern = "RData",
#     recursive = TRUE,
#     full.names = TRUE
# ),
# load,
# .GlobalEnv)

drugbankHmdb <- readRDS(file.path("www", "db", "drugbankHmdb.rds"))
differentialAbundanceSummary <- readRDS(file.path("www", "db", "differentialAbundanceSummary.rds"))

# Convert the fold-change and metdata to log-scale for easier visualization
fc[,1:(dim(fc)[2]-1)] = log2(fc[,1:(dim(fc)[2]-1)])
metdata[,1:(dim(metdata)[2]-2)] = log2(metdata[,1:(dim(metdata)[2]-2)])

# aes_string in ggplotcannot accept names with spaces or colons
#metabs <- paste0("metab", 1:nrow(normal))
#names(metabs) <- rownames(normal)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
        sep="", collapse=" ")
}
