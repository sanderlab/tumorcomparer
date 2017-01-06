# Installation and Use Web App 

    install.packages("devtools")
    
    library(devtools)
    install_bitbucket(repo="cbio_mskcc/tumorcomparer",
        subdir="tumorcomparer",
        auth_user="discoverUser",
        password="discoverUserPassword",
        build_vignette=TRUE,
        dependencies=TRUE,
        args="--no-multiarch")
        
    library(tumorcomparer)
    runShinyApp()

# OTHER
setRepositories(ind=1:6)
options(repos="http://cran.rstudio.com/")

if(!require(devtools)) install.packages("devtools")
if(!require(jsonlite)) install.packages("jsonlite")

library(jsonlite)
library(devtools)

cfg <- '{"repo":"cbio_mskcc/tumorcomparer", "build_vignette":true, "dependencies":true, "auth_user":"discoverUser", "password":"discoverUserPassword", "subdir":"tumorcomparer"}'

cfg <- jsonlite::fromJSON(cfg)
do.call(install_bitbucket, cfg)


git clone http://discoverUser:discoverUserPassword@bitbucket.org/cbio_mskcc/tumorcomparer.git


