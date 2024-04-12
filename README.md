# _seismic_: Single-cell Expression Investigation for Study of Molecular Interactions and Connections
This repository contains code for the R package _seismicGWAS_ a method for calculating cell type-trait associations
given GWAS and single cell RNA-sequencing data.

## Citation

> Disentangling associations between complex traits and cell types with _seismic_.
> Li Q, Dannenfelser R, Roussarie JP, Yao V. BioRxiv. April 2024.

## About

Integrating single-cell RNA sequencing (scRNA-seq) with Genome-Wide Association Studies (GWAS) can help reveal GWAS-associated cell types furthering our understanding of the cell-type-specific biological processes underlying complex traits and disease. In order to rapidly and accurately pinpoint associations, we develop a novel framework, _seismic_, which characterizes cell types using a new specificity score. As part of the _seismic_ framework, the specific genes driving cell type-trait associations can easily be accessed and analyzed, enabling further biological insights. 

![method overview](man/figures/seismic_overview.png)

## Installation
To install the seismic package first clone the seismic repo and then 
use devtools within R to point to seismic and install.

```R
devtools::install(path_to_seismic_folder)
library('seismicGWAS')
```

## Usage
