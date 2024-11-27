# _seismic_: Single-cell Expression Investigation for Study of Molecular Interactions and Connections
This repository contains code for the R package _seismicGWAS_ a method for
calculating cell type-trait associations
given GWAS and single cell RNA-sequencing data.

## Citation

> Disentangling associations between complex traits and cell types with _seismic_.
> Lai Q, Dannenfelser R, Roussarie JP, Yao V. BioRxiv. April 2024.

## About

Integrating single-cell RNA sequencing (scRNA-seq) with Genome-Wide Association
Studies (GWAS) can help reveal GWAS-associated cell types furthering our
understanding of the cell-type-specific biological processes underlying complex
traits and disease. In order to rapidly and accurately pinpoint associations, we
develop a novel framework, _seismic_, which characterizes cell types using a new
specificity score. As part of the _seismic_ framework, the specific genes driving
cell type-trait associations can easily be accessed and analyzed, enabling further
biological insights. The following figure depicts a high level overview of
this process. 

![method overview](man/figures/seismic_overview.png)

## Installation and dependencies
To install the _seismicGWAS_ R package first clone the _seismic_ repo and then 
use devtools within R to point to _seismic_ and install. The R environment version
should be 4.0.0 or higher. 

```bash

```R
devtools::install(path_to_seismic_folder)
library('seismicGWAS')
```

The package require several dependencies which are listed in the `DESCRIPTION` file,
which will be installed automatically.

## Usage
Below we quickly show how to use _seismicGWAS_ to calcuate cell
type-trait associations for the sample data included in the package. 
Full usage instructions, including a walk through of all major functions
can be found in the [vignette](https://github.com/ylaboratory/seismic/tree/master/vignettes/seismicGWAS.md).

```R
# calculate cell type specificity scores using included sample data
# estimate running time: less than 1 minute
tmfacs_sscore <- calc_specificity(tmfacs_sce_small, ct_label_col='cluster_name')

# convert mouse gene identifiers to human ones that match data in GWAS summary data
# from MAGMA
# estimate running time: 3 seconds
tmfacs_sscore_hsa <- translate_gene_ids(tmfacs_sscore, from='mmu_symbol')

# calculate cell type-trait associations for type 2 diabetes
# estimate running time: less than 1 minutes
get_ct_trait_associations(tmfacs_sscore_hsa, t2d_magma)

# find the influential genes for a significant cell type-trait association
# in type 2 diabetes
# estimate running time: 10 seconds
find_inf_genes("Pancreas.beta cell", tmfacs_sscore_hsa, t2d_magma)
```
