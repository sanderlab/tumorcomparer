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