# About

Cell lines derived from human tumors are often used in pre-clinical cancer research, but some cell lines may be too different from tumors to be good models. Genomic and molecular profiles can be used to guide the choice of cell line suitable for particular investigations, but not all features may be equally relevant. We present TumorComparer, a computational method and web service for comparing cell lines and tumors with the flexibility to place a higher weight on functional alterations of interest. In a first pan-cancer application, we compare 260 cell lines and 1914 tumors of six cancer types, using weights emphasizing recurrent genomic alterations. We rank cell lines by their similarity to tumors and identify apparently unsuitable outlier cell lines, including some that are widely used.

# Vignette (Tutorial)

Located in vignettes:

    vignettes/tumor_comparer_example.nb.html

# Installation and Use Web App
    install.packages("devtools")

    library(devtools)
    install_bitbucket(repo="cbio_mskcc/tumorcomparer",
        subdir="tumorcomparer",
        build_vignette=TRUE,
        dependencies=TRUE,
        args="--no-multiarch")

    library(tumorcomparer)
    runShinyApp()

# IGNORE
setRepositories(ind=1:6)
options(repos="http://cran.rstudio.com/")

if(!require(devtools)) install.packages("devtools")
if(!require(jsonlite)) install.packages("jsonlite")

library(jsonlite)
library(devtools)

cfg <- '{"url":"cbio_mskcc/tumorcomparer", "build_vignette":true, "dependencies":true, "subdir":"tumorcomparer"}'

cfg <- jsonlite::fromJSON(cfg)
do.call(install_git, cfg)
