#' Translate gene ids of a given dgeMatrix of seismic specificity scores from a given organism (e.g., mouse)
#' to another (e.g., human) based on orthology.
#'
#' Currently supported organisms: hsa, mmu (human, mouse, respectively)
#' Currently supported IDs: symbol, ensembl, entrez
#'
#' Unmatched rows are dropped from the matrix.
#'
#' @param sscore dgeMatrix of specificity scores (direct output of calc_specificity())
#' @param from The organism and ID type of current identifiers, separated by an underscore (default: mmu_ensembl).
#' @param to The organism and ID type of target identifiers, separated by an underscore (default: hsa_entrez).
#' @param multi_mapping Method to translate metrics when one source ID has several target ID matches. Options: mean (default), sum.
#'
#' @return A dgeMatrix that has row names using target identifiers, with
#' @export
#'
translate_gene_ids <- function(sscore, from = "mmu_ensembl", to = "hsa_entrez", multi_mapping = "mean") {
  # check input type
  if (!class(sscore)[1] %in% c("dgeMatrix")) {
    stop("This function is to be used directly with the direct output of calc_specificity().")
  }

  all_names <- unique(rownames(sscore))

  # currently only supports mmu_hsa_mapping (which is automatically loaded from internal data object sysdata.rda)
  if (!any(all_names %in% mmu_hsa_mapping[[from]])) {
    stop("None of the rownames of sscore could be found in ", from, " identifiers. Please check.")
  }

  # check multi_mapping parameters
  if (!multi_mapping %in% c("mean", "sum")) {
    stop("multi_mapping can only be either 'mean' or 'sum'")
  }

  # filter mapping to only relevant entries
  all_mapping <- unique(mmu_hsa_mapping[, c(from, to), with = F])
  setnames(all_mapping, c("from", "to"))
  all_mapping <- all_mapping[!is.na(from) & !is.na(to) & from %in% all_names]

  # calculate multimapping transformation matrix
  gene_fac_mat <- Matrix::fac2sparse(factor(all_mapping$to, levels = unique(all_mapping$to)))
  if (multi_mapping == "mean") {
    gene_fac_mat <- sweep_sparse(gene_fac_mat, margin = 1, stats = Matrix::rowSums(gene_fac_mat), fun = "/")
  }

  # transform identifiers and return resulting matrix
  return(gene_fac_mat %*% sscore[match(all_mapping$from, rownames(sscore)), ])
}
