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
#' @param gene_mapping_table A data.table, data.frame or tibble with two columns: from and to. 
#' The from column should contain the source gene IDs and the to column should contain the target gene IDs. 
#' The function will use this table to translate the gene IDs. If not provided, the function will use the internal mmu_hsa_mapping data object.
#' @param multi_mapping Method to translate metrics when one source ID has several target ID matches. Options: mean (default), sum.
#'
#' @return A dgeMatrix that has row names using target identifiers, with
#' @export
#'
translate_gene_ids <- function(sscore, from = "mmu_ensembl", to = "hsa_entrez", gene_mapping_table = NULL, multi_mapping = "mean") {
  # check input type
  if (!class(sscore)[1] %in% c("dgeMatrix","dgCMatrix", "matrix")) {
    stop("This function is to be used directly with the direct output of calc_specificity(). The data is not an acceptable matrix type.")
  }

  all_names <- unique(rownames(sscore))
  
  #check gene mapping table
  if (!is.null(gene_mapping_table)) {
    if (!is.data.frame(gene_mapping_table) && !is.data.table(gene_mapping_table) && !is.tibble(gene_mapping_table)) {
      stop("gene_mapping_table must be a data.frame, data.table or tibble")
    }
    gene_mapping_table <- as.data.table(gene_mapping_table)
    #check if the gene_mapping_table contain sufficient information
    if (nrow(gene_mapping_table) < 10000){
      warning("The gene_mapping_table contains less than 10,000 rows.
      There might not be enough information to accurately translate the gene IDs.
      It is recommended to use a more comprehensive gene_mapping table or
      use the default internal mmu_hsa_mapping data object (by not specifying the argument).")
    }
  }else{
    gene_mapping_table <- mmu_hsa_mapping
  }
  
  # currently only supports mmu_hsa_mapping (which is automatically loaded from internal data object sysdata.rda)
  if (!any(all_names %in% gene_mapping_table[[from]])) {
    stop("None of the rownames of sscore could be found in ", from, " identifiers. Please check.")
  }

  #get the unique mapping between from and to
  gene_mapping_table <- gene_mapping_table[!duplicated(gene_mapping_table),]

  #remove rows with NA values
  gene_mapping_table <- gene_mapping_table[!is.na(gene_mapping_table[[from]]) & !is.na(gene_mapping_table[[to]]),]
  
  # check multi_mapping parameters
  if (!multi_mapping %in% c("mean", "sum")) {
    stop("multi_mapping can only be either 'mean' or 'sum'")
  }
 
  # filter mapping to only relevant entries
  all_mapping <- unique(gene_mapping_table[, c(from, to), with = F])
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
