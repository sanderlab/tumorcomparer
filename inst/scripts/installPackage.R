setRepositories(ind=1:6)
options(repos="http://cran.rstudio.com/")

if(!require(devtools)) { install.packages("devtools") }

library(devtools);

install_bitbucket("cbio_mskcc/tumorcomparer",
                  build_vignette=FALSE,
                  dependencies=TRUE,
                  args="--no-multiarch",
                  subdir="tumorcomparer")
