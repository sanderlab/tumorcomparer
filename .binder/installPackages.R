#' Install/update necessary packages from CRAN, Bioconductor, GitHub, or local sources
#'
#' @param file a file with packages; overrides packages parameter
#' @param packages a vector of strings with names of packages from CRAN, Bioconductor, GitHub
#' @param updatePackages whether to update existing packages (Default: FALSE)
#' @param dryRun whether to test for missing packages (Default: FALSE)
#'
#' @example 
#' \dontrun {
#' source("https://gist.githubusercontent.com/cannin/6b8c68e7db19c4902459/raw/installPackages.R")
#' installPackages("r-requirements.dcf")
#' }
installPackages <- function(file=NULL, packages=NULL,
                            type=getOption("pkgType"), 
                            repos="http://cran.rstudio.com/",
                            updatePackages=FALSE,
                            buildGitVignettes=FALSE,
                            dryRun=FALSE) {
  setRepositories(ind=1:6)
  options(repos=repos, unzip="internal") # unzip needed for Ubuntu
  
  # Install required packages 
  if (!requireNamespace("remotes", quietly = TRUE)) { install.packages("remotes") }
  if (!requireNamespace("git2r", quietly = TRUE)) { install.packages("git2r") }
  
  if(is.null(packages) && is.null(file)) {
    stop("ERROR: Either packages or file must be set.")
  }
  
  # Read requirements DCF 
  if(!is.null(file)) {
    dcf <- read.dcf(file)
    
    tryCatch({
      rVersion <- numeric_version(dcf[1, "r-version"])
      rVersion <- numeric_version(paste0(rVersion$major, rVersion$minor))
      
      rVersionSystem <- getRversion()
      rVersionSystem <- numeric_version(paste0(rVersionSystem$major, rVersionSystem$minor))
      
      # Check R Version 
      stopifnot(rVersionSystem == rVersion)
    }, error = function(e) {
      cat("WARNING: Optional r-version parameter not set.\n")
    })
    
    tryCatch({
      mranDate <- dcf[1, "mran-date"]
      mranRepos <- paste('https://mran.microsoft.com/snapshot/',  mranDate, '/', sep="")
      
      # Set MRAN Repos 
      options(repos=mranRepos, unzip="internal")      
    }, error = function(e) {
      cat("WARNING: Optional mran-date parameter not set.\n")
    })

    tryCatch({
      updatePackages <- dcf[1, "update-packages"]
      updatePackages <- as.logical(updatePackages)      
    }, error = function(e) {
      cat("WARNING: Optional update-packages parameter not set.\n")
    })

    depends <- dcf[1, "depends"]
    packages <- strsplit(depends, "\n")[[1]]
  }
  
  # Install remotes 
  if("remotes" %in% rownames(installed.packages())) {
    require(remotes)      
  } else {
    install.packages("remotes")  
    require(remotes)    
  }
  
  # Install packages 
  if(!dryRun) {
    for(package in packages) {
      cat("Processing: ", package, "\n")
      packageName <- package 
      
      # Get just the package name 
      if(grepl("::", package)) {
        packageName <- strsplit(package, "::")[[1]][2]
      }
      packageName <- strsplit(packageName, .Platform$file.sep)[[1]]
      packageName <- packageName[length(packageName)]
      
      if(!(packageName %in% rownames(installed.packages()))) {
        tryCatch({
          if(package %in% rownames(available.packages())) {
            install_cran(package, type=type, upgrade=updatePackages)
          } else if(grepl("^github::", package)) {
            tmpPkg <- strsplit(package, "::")[[1]][2]
            install_github(tmpPkg, upgrade=updatePackages, build_vignettes=buildGitVignettes)
          } else if(grepl("^bitbucket::", package)) {
            tmpPkg <- strsplit(package, "::")[[1]][2]
            install_bitbucket(tmpPkg, upgrade=updatePackages, build_vignettes=buildGitVignettes)
          } else if(grepl("^local::", package)) {
            tmpPkg <- strsplit(package, "::")[[1]][2]
            install_local(tmpPkg, repos = NULL, type="source")
          } else if(grepl("^bioc::", package)) {
            tmpPkg <- strsplit(package, "::")[[1]][2]
            install_bioc(tmpPkg, upgrade=updatePackages)
          } else {
            cat("WARNING: No remote:: specified. Trying package: ", package, "with Bioconductor.\n")
            if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
            BiocManager::install(package, update=FALSE, ask=FALSE)
          }
          
          # Try to test if rJava worked
          # NOTE: Requires additional steps: https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite
          if(packageName == "rJava") {
            library(rJava) 
            .jinit()
          }
        }, error = function(e) {
          cat("ERROR: Package: ", package, ". Message: ", message(e), "\n")
        })
      } else {
        cat("Already installed: ", package, "\n")
      }
    }    
  } 

  packageNames <- sapply(packages, function(package) {
    if(grepl("::", package)) {
      package <- strsplit(package, "::")[[1]][2]
    }
    
    package <- strsplit(package, .Platform$file.sep)[[1]]

    packageName <- package[length(package)]
  }, USE.NAMES=FALSE)
  
  idx <- which(!(packageNames %in% rownames(installed.packages())))
  cat("Missing packages: ", paste(packageNames[idx], collapse=", "), "\n")
}
