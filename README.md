# About

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
