# Table of Contents 

- [About TumorComparer](#about-tumorcomparer)
- [Analysis](#analysis)
  * [Pre-Computed Analysis](#pre-computed-analysis)
  * [User Data Analysis](#user-data-analysis)
    + [Parameters](#parameters)
    + [Data Formats](#data-formats)
      - [Expression Data](#expression-data)
      - [Mutation Data](#mutation-data)
      - [Copy Number Data](#copy-number-data)
  * [Example Datasets](#example-datasets)
- [Additional Features](#additional-features)
- [Feedback](#feedback)

# About TumorComparer
Built by Rileen Sinha, Augustin Luna, Nikolaus Schultz, and Chris Sander using R Shiny.

FIXME: Cell lines derived from human tumors are often used in pre-clinical cancer research, but some cell lines may be too different from tumors to be good models. Genomic and molecular profiles can be used to guide the choice of cell line suitable for particular investigations, but not all features may be equally relevant. We present TumorComparer, a computational method and web service for comparing cell lines and tumors with the flexibility to place a higher weight on functional alterations of interest. In a first pan-cancer application, we compare 260 cell lines and 1914 tumors of six cancer types, using weights emphasizing recurrent genomic alterations. We rank cell lines by their similarity to tumors and identify apparently unsuitable outlier cell lines, including some that are widely used.

# Analysis

The "Pre-Computed" tab allows users to explore the results of our systematic analysis available in the publication, while the "User Analysis" section provides users the option of running TumorComparer with their own data.

## Pre-Computed Analysis

All cancer types from the publication are available. Balloon plots show the most relevant cell lines and a corresponding searchable table.

## User Data Analysis

**FIXME**

### Parameters 

* Default (Background) Weight: 0.01 (DEFAULT); all other genes

Users are should review the publication for more information about these parameters. The default weights were used as part of the systematic analysis. 

Users should rely on their understanding of the problem they are trying to address to select weights. Users can use 1 and 0 as a starting point for known cancer gene and default weights if the relative difference in the magnitude in importance between the sets of genes is unknown.

### Data Formats 

Below is sample data for 

#### Expression Data

|  |RCM.1|SW837|SW1463|SNU.61|CaR.1|
|--------|-----|-----|------|------|-----|
|AREG    |11.837|9.195|9.018 |9.939 |10.44|
|MLXIPL  |9.513|8.823|11.749|9.008 |8.029|
|MUC16   |3.062|7.281|6.235 |2.774 |9.931|
|ADAMTS15|6.177|5.697|5.272 |6.198 |6.004|
|HBB     |7.088|2.609|1.093 |1.379 |6.138|
|COL1A1  |6.967|6.647|5.955 |5.993 |6.891|


#### Mutation Data

|  |RCM.1|SW837|SW1463|SNU.61|CaR.1|
|--------|-----|-----|------|------|-----|
|ANKRD30A|0    |0    |0     |0     |0    |
|APC     |1    |1    |1     |1     |0    |
|APOB    |0    |0    |0     |0     |0    |
|ARHGAP35|0    |0    |0     |0     |0    |
|ARID1A  |0    |0    |0     |0     |0    |
|ARID1B  |0    |0    |0     |0     |0    |


#### Copy Number Data

Below is an example of GISTIC discretized 

|  |RCM.1|SW837|SW1463|SNU.61|CaR.1|
|--------|-----|-----|------|------|-----|
|PFN1P2  |0    |0    |0     |0     |0.667|
|PDE4DIP |0    |0    |0     |0     |0.667|
|SEC22B  |0    |0    |0     |0     |0.667|
|NOTCH2NL|0    |0    |0     |0     |0.667|
|NBPF10  |0    |0    |0     |0     |0.667|
|HFE2    |0    |0    |0     |0     |0.667|

## Example Datasets 

A number of example datasets are available:

* Rectal Adenocarcinoma (READ) TCGA/CCLP: Small dataset, [LINK]()
* Ovarian (OVCA) TCGA/CCLP: Larger dataset, [LINK]()
* Melanoma TCGA/CCLP, Custom Gene Lists: Provides users an example for creating customized lists based on pathway genes taken from the [TCGA PanCan Pathways Marker Paper](https://pubmed.ncbi.nlm.nih.gov/29625050/).
 * RTK-RAS Pathway: [LINK]()
 * WNT Pathway: [LINK]()
 
# Additional Features 

The TumorComparer software package contains a number of additional parameters that might be of interest to users. 

# Feedback
We appreciate any feedback/suggestions you may have; please feedback to publication corresponding authors.
