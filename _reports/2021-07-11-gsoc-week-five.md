---
layout: post
title:  "Week five | TumorComparer on Gene Subsets"
tags: [gsoc, weekly report, coding period ]
author: Siras Hakobyan
---

This week I was calculating comparisons for 10 specified genesets. During the calculation, we detected some technical issues and exceptions which were not covered by package functions. With the help of my mentor, I have updated package functions to handle exceptions and missing values. 

## Work Progress

1. **TumorComparer on Gene Subsets**   
    Status։ **Done**    
    **Issue:** [TumorComparer on Gene Subsets](https://github.com/sanderlab/tumorcomparer/issues/11)

2. **Reproducibility of the results**   
    Status։ **Done**    
    **Issue:** [Testing Geneset Specific Function](https://github.com/sanderlab/tumorcomparer/issues/13)

2. **Missing values and exceptions**   
    Status։ **Done**    
    **Issue:** [Geneset filtering issue](https://github.com/sanderlab/tumorcomparer/issues/15)

3. **Missing values and exceptions**   
    Status։ **Done**    
    **Issue:** [Datatypes with zero values in precomputed comparisons](https://github.com/sanderlab/tumorcomparer/issues/16)

4. **Missing values and exceptions**   
    Status։ **Done**    
    **Issue:** [Add Test that Triggers to Few Genes in Comparison](https://github.com/sanderlab/tumorcomparer/issues/17)

5. **Missing values and exceptions**   
    Status։ **Done**    
    **Issue:** [Comparison function error with LIHC cancer type in TGF-Beta pathway](https://github.com/sanderlab/tumorcomparer/issues/18)


## Conclusion  

I have finished precalculation of comparisons for 10 genesets in 23 cancer types, updated package functions to handle missing and zero values. On top of that, I have updated the generalized comparison function to perform comparison for specified gene list.
