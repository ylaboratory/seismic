#' Calculate influential genes for a given trait and cell type using DFBETAS.
#'
#' @param ct A character string containing a valid cell type name in sscore.
#' @param sscore A dgeMatrix of seismic specificity scores where
#' each column is a cell type and row names are gene identifiers.
#' (Note: the identifiers used should match those used in the MAGMA input)
#' @param magma A data.frame or file path to MAGMA output for a particular GWAS
#' with at least 2 columns: gene identifiers and z-scores.
#' @param magma_gene_col A character string corresponding to the column name
#' in the MAGMA data containing gene identifiers. Defaults to 'GENE'
#' from the MAGMA output.
#' @param magma_z_col A character string corresponding to the column name
#' in the MAGMA data containing z-scores. Defaults to 'ZSTAT'
#' from the MAGMA output.
#'
#' @return A data.frame containing genes, seismic specificity scores, magma
#' z-stats, the dfbetas influential gene score, and a Boolean value indicating
#' if the gene is influential.
#' @export
find_inf_genes <- function(ct, sscore, magma,
                           magma_gene_col = "GENE", magma_z_col = "ZSTAT") {
  dfbetas <- is_influential <- NULL # due to non-standard evaluation notes in R CMD check

  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)

  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)

  # check that cell type name exists
  if (!(ct %in% colnames(sscore))) {
    stop("Invalid cell type name. Please choose a cell type
         present in the given specificity data.")
  }

  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- sscore[, c("gene", ct), with = F]
  names(sscore) <- c("gene", "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$specificity <- as.numeric(sscore$specificity)

  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)

  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]
  lm_out <- stats::lm(zstat ~ specificity, data = dt)
  lm_summ <- summary(lm_out)$coefficients
  lm_pval <- if (lm_summ[2, 1] > 0) lm_summ[2, 4] / 2 else (1 - lm_summ[2, 4] / 2)
  if (lm_pval > 0.05) {
    warning("Cell type ", ct, " does not seem to have a significant association with trait
         (estimated p-value ", round(lm_pval, 4), " based on a 1-sided hypothesis test,
         prior to multiple hypothesis test correction, so influential gene analysis may not be as relevant")
  }

  dt[, dfbetas := stats::dfbetas(lm_out)[, 2]]
  thresh <- 2 / sqrt(nrow(dt))
  # only considering positive dfbetas values (since our analyses is based off of
  # the 1-sided test of positive relationships between specificity + risk)
  dt[, is_influential := ifelse(dfbetas > thresh, T, F)]
  dt <- dt[order(-dfbetas)]

  return(dt)
}
