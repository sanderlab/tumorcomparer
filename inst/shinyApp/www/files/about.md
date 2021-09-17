# Table of Contents 

- [About TumorComparer](#about-tumorcomparer)
- [Analysis](#analysis)
  * [Pre-Computed Analysis](#pre-computed-analysis)
  * [User Data Analysis](#user-data-analysis)
    + [Parameters](#parameters)
    + [Data Formats](#data-formats)
      - [Input files](#input-files)
      - [Example Data](#example-data)
        * [Notes](#notes)
      - [Expression Data](#expression-data)
      - [Mutation Data](#mutation-data)
      - [Copy Number Data](#copy-number-data)
  * [Example Datasets](#example-datasets)
- [Additional Features](#additional-features)
- [Feedback](#feedback)

# About TumorComparer
## Team 
Built by Rileen Sinha, Augustin Luna, Nikolaus Schultz, and Chris Sander. 

# Abstract
Cancer is a genetic disease, typically marked by widespread somatic alterations (e.g., mutations, copy-number alterations, and gene expression changes). However, not all changes are functionally important—few genes can promote oncogenesis (also termed “cancer drivers”), whereas other altered genes have little effect on the phenotype (termed “passengers”). Furthermore, many research questions focus on particular genes and their activity (e.g., specific signaling pathways, drug targets, etc.). This motivates the need for a flexible method of comparing tumors with potential cell line models by using researcher-selected properties. We present TumorComparer, a computational comparison method based on weighted features to allow expert- and knowledge-driven comparison of tumors and experimental models, such as cell lines or organoids. We apply TumorComparer to the comparison of ∼8,000 tumors and ∼600 cell lines across 24 cancer types as an initial application to provide a general, pan-cancer resource based on knowledge of oncogenic alterations gained from The Cancer Genome Atlas program (TCGA). TumorComparer is a generally applicable method suitable for pre-clinical cancer research and personalized medicine applications where sets of samples need to be assessed for similarity.

## Publication
https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(21)00084-9

# Analysis

The "Pre-Computed" tab allows users to explore the results of our systematic analysis available in the publication, while the "User Analysis" section provides users the option of running TumorComparer with their own data.

## Pre-Computed Analysis

All cancer types from the publication are available. Balloon plots show the most relevant cell lines and a corresponding searchable table.

## User Data Analysis

This tab allows users to run TumorComparer on their own data. Users will upload a zip file containing the input files with a pre-specified naming convention. 

NOTE: For large computational pipelines or large amounts of data, it is suggested that users run TumorComparer locally. 

### Parameters 

* Default (Background) Weight: 0.01 (DEFAULT); all other genes

Users are should review the publication for more information about these parameters. The default weights were used as part of the systematic analysis. 

Users should rely on their understanding of the problem they are trying to address to select weights. Users can use 1 and 0 as a starting point for known cancer gene and default weights if the relative difference in the magnitude in importance between the sets of genes is unknown.

### Data Formats

#### Input files 

NOTE: You do NOT have to provide input files for all data types (expression: exp, mutation: mut, copy number: cna). Make sure either ALL files for a data type are included OR NONE of them are included, based on your available data.

* tumor_mut.txt: a file with binary mutation data for tumors 
* tumor_cna.txt: a file with GISTIC data for tumors; this can be 5-values (-2, -1, 0, 1, 2) or continuous
* tumor_exp.txt: a file with gene expression data for tumors
* cell_line_mut.txt: See corresponding tumor file 
* cell_line_exp.txt: See corresponding tumor file 
* cell_line_cna.txt: See corresponding tumor file
* genes_and_weights_mut.txt: a file with weights for cancer-specific set of recurrently mutated genes. A tab-delimited file - the first column has the gene names, and the second column specifies the weights.
* genes_and_weights_cna.txt: See genes_and_weights_mut.txt
* genes_and_weights_exp.txt: See genes_and_weights_mut.txt
* default_weights_for_known_cancer_genes_mut.txt: a file with weights for genes known to be recurrently altered/mutated in cancer (e.g. recurrently mutated genes in TCGA pan-cancer analyses). A two-column tab-delimited file - the first column has the gene names and the second column specifies the weights.
* default_weights_for_known_cancer_genes_exp.txt: See default_weights_for_known_cancer_genes_mut.txt
* default_weights_for_known_cancer_genes_cna.txt: See default_weights_for_known_cancer_genes_mut.txt

#### Example Data

Below is example data for the various supported data types. 

##### Notes

* Sample names should match within the tumor data and separately, within the cell line data
* Analysis is platform agnostic (see publication; i.e., microarray or RNA-seq can be used)

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

* Rectal Adenocarcinoma (READ) TCGA/CCLP: Small dataset, [Link](./read_data_for_running_tc.zip)

# Additional Features 

The TumorComparer software package contains a number of additional parameters that might be of interest to users. 

# Feedback

We appreciate any feedback/suggestions you may have; please send feedback to publication corresponding authors.
