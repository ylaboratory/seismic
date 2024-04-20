#' Small Subset of Tabula Muris FACS scRNA-seq data
#'
#' A subset of 4000 cells from the Tabula Muris FACS scRNA-seq data
#' that has already been normalized using the `scran` package, providing
#' the log-normalized counts. Cluster labels are under a labeled `cluster_name` column.
#'
#' @format ## `tmfacs_sce_small`
#' A SingleCellExperiment with 23,341 rows (mouse gene symbols) and 4000 columns (cells),
#' with relevant items:
#' \describe{
#'   \item{assays$logcounts}{log normalized counts calculated from scran}
#'   \item{cluster_name}{tissue and cell type (separated by .)}
#' }
#' @source <https://tabula-muris.ds.czbiohub.org>
"tmfacs_sce_small"

#' MAGMA genes.out output for rheumatoid arthritis
#'
#' A data.frame of the `genes.out` output from MAGMA calculated using
#' rheumatoid arthritis summary statistics.
#'
#' @format ## `ra_magma`
#' A data.frame with 17,460 rows and 9 columns:
#' \describe{
#'   \item{GENE}{entrez gene ID}
#'   \item{CHR}{the chromosome the gene is on}
#'   \item{START}{start boundary of the gene on the chromosome (including any provided window during annotation)}
#'   \item{END}{end boundary of the gene on the chromosome (including any provided window during annotation)}
#'   \item{NSNPS}{number of SNPS annotated to gene (not excluded based on internal SNP QC)}
#'   \item{NPARAM}{number of relevant parameters used in the MAGMA model}
#'   \item{N}{sample size used in analyzing the gene}
#'   \item{ZSTAT}{z-value of gene based on permuted p-value (what we also term MAGMA z-score)}
#'   \item{P}{gene p-value}
#' }
#' @source <https://www.nature.com/articles/nature12873>
"ra_magma"

#' MAGMA genes.out output for type 2 diabetes
#'
#' A data.frame of the `genes.out` output from MAGMA calculated using
#' type 2 diabetes summary statistics.
#'
#' @format ## `t2d_magma`
#' A data.frame with 17,460 rows and 9 columns:
#' \describe{
#'   \item{GENE}{entrez gene ID}
#'   \item{CHR}{the chromosome the gene is on}
#'   \item{START}{start boundary of the gene on the chromosome (including any provided window during annotation)}
#'   \item{END}{end boundary of the gene on the chromosome (including any provided window during annotation)}
#'   \item{NSNPS}{number of SNPS annotated to gene (not excluded based on internal SNP QC)}
#'   \item{NPARAM}{number of relevant parameters used in the MAGMA model}
#'   \item{N}{sample size used in analyzing the gene}
#'   \item{ZSTAT}{z-value of gene based on permuted p-value (what we also term MAGMA z-score)}
#'   \item{P}{gene p-value}
#' }
#' @source <https://www.nature.com/articles/s41467-018-04951-w>
"t2d_magma"
