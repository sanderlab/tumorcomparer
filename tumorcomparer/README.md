# Installation and Usage 

    install.packages("devtools")
    
    library(devtools)
    install_bitbucket("cbio_mskcc/pancanmet",
        auth_user="discoverUser",
        password="discoverUserPassword",
        build_vignette=TRUE,
        dependencies=TRUE,
        args="--no-multiarch")
        
    library(pancanmet)
    runShinyApp()
    
#' @param composite_mut a composite (both tumor and cell line information) matrix 
#' with only mutation information
#' @param composite_cna a composite (both tumor and cell line information) matrix 
#' with only copy number information