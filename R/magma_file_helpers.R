#' Check magma object / file is usable and returns it as a data.table
#'
#' @param magma input magma object / file
#' @param magma_gene_col character string corresponding to the column name
#' in the MAGMA data containing gene identifiers
#' @param magma_z_col character string corresponding to the column name
#' in the MAGMA data containing z-scores
#'
#' @return magma file contents as a data.table
load_magma_dt <- function(magma, magma_gene_col, magma_z_col) {
  if (is.character(magma)) {
    if (!file.exists(magma)) {
      stop("input magma file ", magma, " does not exist")
    }
    magma <- fread(magma)
  }

  if (is.data.frame(magma)) { # will return T for data.frame or data.table etc
    magma <- as.data.table(magma)
  }

  if (!(magma_gene_col) %in% names(magma)) {
    stop("magma_gene_col ", magma_gene_col, " cannot be found in input magma file")
  }

  if (!(magma_z_col) %in% names(magma)) {
    stop("magma_z_col ", magma_z_col, " cannot be found in input magma file")
  }

  return(magma)
}


#' Check for overlapping genes between specificity scores and trait gene
#' summaries from MAGMA.
#'
#' @param sscore A dgeMatrix of seismic specificity scores where
#' each column is a cell type and row names are gene identifiers.
#' (Note: the identifiers used should match those used in the MAGMA input)
#' @param magma A data.frame or file path to MAGMA output for a particular GWAS
#' with at least 2 columns: gene identifiers and z-scores.
#' @param magma_gene_col A character string corresponding to the column name
#' in the MAGMA data containing gene identifiers. Defaults to 'GENE'
#' from the MAGMA output.
#' @param warn_thresh A floating point number corresponding to the minimum
#' percentage of overlapping genes before the warning error is displayed.
#' Defaults to 80% overlap.
#' @param min_intersect A floating point number representing the minimum number
#' of overlapping genes needed for further execution. Defaults to 1.
#'
check_overlap <- function(sscore, magma, magma_gene_col, warn_thresh = 0.8,
                          min_intersect = 1) {
  tot_genes <- length(rownames(sscore))
  int_genes <- length(intersect(
    rownames(sscore),
    as.character(magma[[magma_gene_col]])
  ))

  if ((int_genes / tot_genes) < warn_thresh) {
    warning("Only ", (int_genes / tot_genes) * 100, "% (", int_genes, " genes)
          of genes map between seismic specificity scores and the MAGMA
          input. If this is unexpected, please check the identifiers between
          these files or change the magma_gene_col and try again.")
  }

  if (int_genes < min_intersect) {
    stop("Only ", min_intersect, " overlapping gene(s) between seismic specificity scores
         and the MAGMA input. This is below the minimum threshold needed to run.
         Please check the identifiers between
         these files or change the magma_gene_col and try again.")
  }
}
