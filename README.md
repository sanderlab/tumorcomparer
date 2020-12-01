# Abstract

Cell lines derived from human tumors are often used in pre-clinical cancer research, but some cell lines may be too different from tumors to be good models. Genomic and molecular profiles can be used to guide the choice of cell lines suitable for particular investigations, but not all features may be equally relevant. We present TumorComparer, a computational method for comparing cellular profiles  with the flexibility to place a higher weight on functional features of interest. In a first pan-cancer application, we compare 600 cell lines and 8323 tumors of 26 cancer types, using weights to emphasize recurrent genomic alterations or expression dysregulation. We characterize the similarity of cell lines and tumors within and across cancers, identifying outlier and mislabelled cell lines as well as good matches, and identify cancers with an unusually high number of good or poor representative cell lines. The weighted similarity method in the future may be useful to assess genomic-molecular patient profiles for stratification in clinical trials and personalized choice of therapy.

# Installation

```
install.packages("devtools")

library(devtools)
install_bitbucket(repo="cbio_mskcc/tumorcomparer",
    build_vignette=TRUE,
    dependencies=TRUE,
    args="--no-multiarch")
```

# Vignette (Tutorial)

A pre-built vignette is committed at 'vignettes/tumor_comparer_example.nb.html'

# Run Web Application Locally 

```
library(tumorcomparer)
runShinyApp()
```

