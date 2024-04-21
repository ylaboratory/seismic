#' Calculate cell type-trait associations
#'
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
#' @return A data.frame of trait associations for each cell type.
#' @export
get_ct_trait_associations <- function(sscore, magma, magma_gene_col = "GENE",
                                      magma_z_col = "ZSTAT") {
  pvalue <- FDR <- NULL # due to non-standard evaluation notes in R CMD check

  # if given file path to magma, check that it is loadable
  magma <- load_magma_dt(magma, magma_gene_col, magma_z_col)

  # check that there is overlap enough to run
  check_overlap(sscore, magma, magma_gene_col)

  # clean data formatting
  sscore <- as.data.table(as.matrix(sscore), keep.rownames = T)
  setnames(sscore, "rn", "gene")
  sscore <- melt(sscore, id.vars = "gene", variable.name = "cell_type", value.name = "specificity")
  sscore$gene <- as.character(sscore$gene)
  sscore$cell_type <- as.character(sscore$cell_type)
  sscore$specificity <- as.numeric(sscore$specificity)

  magma <- magma[, c(magma_gene_col, magma_z_col), with = F]
  names(magma) <- c("gene", "zstat")
  magma$gene <- as.character(magma$gene)

  # get all cell type names
  cts <- unique(sscore$cell_type)

  # combine the magma annots with the specificities
  dt <- merge(sscore, magma, by = "gene")
  dt <- dt[stats::complete.cases(dt)]

  # calculate the association for each cell type
  res <- rbindlist(lapply(cts, function(ct) {
    slm <- speedglm::speedlm(dt[cell_type == ct]$zstat ~ dt[cell_type == ct]$specificity) # fast lm
    slm_summ <- summary(slm)$coefficients

    # extract the 1-sided hypothesis test p-value
    pval <- if (slm_summ[2, 1] > 0) slm_summ[2, 4] / 2 else (1 - slm_summ[2, 4] / 2)
    return(data.table(cell_type = ct, pvalue = pval))
  }))

  res[, FDR := stats::p.adjust(pvalue, method = "fdr")]
  res <- res[order(pvalue, FDR)]

  return(res)
}
