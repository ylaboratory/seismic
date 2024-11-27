#' Compute seismic specificity score for each gene and cell type for a list of SingleCellExperiment objects.
#'
#' @param sce_list A list of multiple SingleCellExperiment objects. This object needs to include the
#' assay_name specified and a ct_label_col column corresponding to cell type labels
#' for the granularity of interest. Row names are used as the gene name identifiers.
#' @param assay_name An assay in sce, default 'logcounts' (calculated using scran).
#' @param ct_label_col A column name in colData of sce (information for the cells in sce)
#' that provides labels for each cell, indicating membership in each cell type
#' at the granularity of interest.
#' @param min_uniq_ct The minimum number of unique cell types needed in the data to
#' proceed with specificity score calculations. Default: 2.
#' @param min_ct_size Filter for cell types to be included in the specificity score calculation
#' based on the number of cells that are labeled as that cell type. Default: 20.
#' @param min_cells_gene_exp Filter for genes to be included in the specificity score calculation
#' based on the number of cells where the gene has non-zero expression. Default: 10.
#' @param min_avg_exp_ct Filter for genes to be included in the specificity score calculation
#' based on the mean of cell-type-averaged expressions (aka mean of mean cell type expression). Default: 0.1
#' @param out_group_mat A matrix indicating the outgroup information -> 
#'
#' @return A dgeMatrix of specificity scores (genes as rows, cell types as columns).
#' @export
calc_specificity_multi_obj <- function(sce_list, assay_name = "logcounts",
                             ct_label_col = "idents", 
                             out_group_mat = NULL,
                             min_uniq_ct = 2,
                             min_ct_size = 20, 
                             min_cells_gene_exp = 10,
                             min_avg_exp_ct = 0.1) {
  ct <- N <- nz.count <- ave_exp_ct <- NULL # due to non-standard evaluation notes in R CMD check

  # check if input is a list of SingleCellExperiment objects
  if (!is.list(sce_list) || length(sce_list) == 0 ||
        !all(sapply(sce_list, function(x) inherits(x, "SingleCellExperiment")))) {
    stop("Input must be a list of SingleCellExperiment objects.")
  }
  
  # check if the assay_name and ct_label_col are present in all SingleCellExperiment objects
  if (!all(sapply(sce_list, function(x) assay_name %in% SummarizedExperiment::assayNames(x)))) {
    stop("Assay '", assay_name, "' does not exist in at least one element
         of the SingleCellExperiment object list. Please choose a valid
         assay in the SingleCellExperiment.")
  }
  

  if (!all(sapply(sce_list, function(x) ct_label_col %in% colnames(SummarizedExperiment::colData(x))))){
    stop("ct_label_col '", ct_label_col, "'does not exist in at least one element
         of the SingleCellExperiment object list. Please choose a valid
         column of cell type labels in the SingleCellExperiment.")
  }

  #check if they have the same rownames
  if (!all(sapply(sce_list, function(x) identical(SummarizedExperiment::rownames(sce_list[[1]]), SummarizedExperiment::rownames(x))))){
    stop("The SingleCellExperiment objects in the list do not have the same rownames.")
  }

  # extract the log normalized counts (dgCMatrix)
  data_mat_list <- lapply(sce_list, function(x) SummarizedExperiment::assay(x, assay_name))


  #check if the matrices have rownames
  if (any(sapply(data_mat_list, function(x) any(is.null(rownames(x)))))) {
    stop("Some data matrices of the assay in the SingleCellExperiment object list do not have rownames.")
  }

  # check that the data matrices have the same rownames
  if (!all(sapply(data_mat_list, function(x) identical(rownames(data_mat_list[[1]]), rownames(x))))){
    stop("The data matrices in the SingleCellExperiment objects do not have the same rownames.")
  }
  
  
  # extract the cell metadata
  cell_meta_list <- lapply(sce_list, function(x) SummarizedExperiment::colData(x)) 

  # make sure the cell name exists in the cell metadata
  init_idx = 1
  for (i in 1:length(cell_meta_list)){
    if (any(is.null(rownames(cell_meta_list[[i]]))) || any(is.null(colnames(data_mat_list[[i]])))){
      rownames(cell_meta_list[[i]]) <- colnames(data_mat_list[[i]]) <- paste0("cell.",init_idx: (ncol(data_mat_list[[i]] ) + init_idx - 1)) 
      init_idx = init_idx + ncol(data_mat_list[[i]])
    }
  }

  #check if cell indeces overlap
  if (length(unique(unlist(lapply(cell_meta_list, rownames)))) != sum(unlist(lapply(data_mat_list, ncol)))){
    stop("The cell names in the cell metadata are not unique.")
  }
  
  # extract cell type grouping
  ct_groups_list <- lapply(cell_meta_list, function(x) data.table(
    cell = rownames(x),
    ct = x[[ct_label_col]], key = "ct"
  ))

  # check that there are at least a few different cell types
  ct_groups_n_list <- lapply(ct_groups_list, function(x) x[, .N, by = ct])

  #summarised N
  ct_groups_n <- rbindlist(ct_groups_n_list)[,.(N = sum(N)), by = ct]

  if (nrow(ct_groups_n) < min_uniq_ct) {
    stop("There are fewer than ", min_uniq_ct, " in the SingleCellExperiment.
         Decrease the min_uniq_ct threshold or select cell type labels with
         more unique cell types.")
  }

  # filter out cell types that do not have minimum # of cells
  ct_groups_n <- ct_groups_n[N >= min_ct_size]
  ct_groups_n_list <- lapply(ct_groups_n_list, function(x) x[ct %in% ct_groups_n$ct])
  ct_groups_list <- lapply(ct_groups_list, function(x) x[ct %in% ct_groups_n$ct])
  data_mat_list <- lapply(seq_along(data_mat_list), function(i) data_mat_list[[i]][, ct_groups_list[[i]]$cell])
  print("line 111")
  if (nrow(ct_groups_n) < min_uniq_ct) {
    stop("There are fewer than ", min_uniq_ct, " in the SingleCellExperiment after
          filtering out cell types with fewer than ", min_ct_size, " cells.
         Decrease the min_uniq_ct or min_ct_size thresholds or select cell type labels with
         more unique cell types / larger numbers of cells.")
  }

  # filter out genes that do not have minimal cell coverage
  stats.dt <- lapply(data_mat_list, function(x) as.data.table(Matrix::rowSums(x != 0), keep.rownames = T) %>%
    magrittr::set_colnames(c("gene", "nz.count"))) %>%
    rbindlist() %>%
    .[, .(nz.count = sum(nz.count)), by = gene]
  
  #stats.dt <- stats.dt[nz.count >= min_cells_gene_exp]

  if (nrow(stats.dt) < 1000) {
    warning("There are fewer than 1000 genes in the SingleCellExperiment that
            are expressed in at least ", min_cells_gene_exp, " cells.
            Consider relaxing the threshold or double check the input file.")
  }

  data_mat_list <- lapply(data_mat_list, function(x) x[stats.dt$gene, ])

  # calculate mean gene expression per cell type
  factor_mat_list <- lapply(ct_groups_list, function(x) Matrix::fac2sparse(factor(x$ct, levels = unique(ct_groups_n$ct)), drop.unused.levels = FALSE))
  
  sum_mat <- lapply(seq_along(data_mat_list), function(i) Matrix::t(data_mat_list[[i]] %*% Matrix::t(factor_mat_list[[i]]))) %>% 
    Reduce(`+`, .) 

  mean_mat <- sum_mat %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat_list[[1]])) %>%
    magrittr::set_rownames(unique(ct_groups_n$ct))
  
  print("line 146")
  
  # filter out genes whose cell-type averaged expression does not exceed a baseline level
  stats.dt$ave_exp_ct <- Matrix::colMeans(mean_mat)
  #stats.dt <- stats.dt[ave_exp_ct >= min_avg_exp_ct]

  if (nrow(stats.dt) < 1000) {
    warning("There are fewer than 1000 genes in the SingleCellExperiment that
            do not have at least mean expression ", min_avg_exp_ct, " across cell types.
            Consider relaxing the threshold or double check the input file.")
  }

  data_mat_list <- lapply(data_mat_list, function(x) x[stats.dt$gene, ])
  sum_mat <- sum_mat[, stats.dt$gene]
  mean_mat <- mean_mat[, stats.dt$gene]
  print("line 161")
  # calculate variance of gene expression per cell type
var_mat <- lapply(seq_along(data_mat_list), function(i) (Matrix::t(data_mat_list[[i]]^2 %*% Matrix::t(factor_mat_list[[i]])))) %>% # calculate the sum of squares divided by N-1
    Reduce(`+`, .) %>%
    {. - 2 * mean_mat * sum_mat} %>%
    {. + sweep_sparse(mean_mat^2, margin = 1, stats = ct_groups_n$N, fun = "*")} %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat_list[[1]])) %>%
    magrittr::set_rownames(unique(ct_groups_n$ct))
  print("line 170")
  # calculate ratio of cells a gene has non-zero expression in per cell type
  ratio_mat <- lapply(seq_along(data_mat_list), function(i) Matrix::t((data_mat_list[[i]] > 0) %*% Matrix::t(factor_mat_list[[i]]))) %>%
    Reduce(`+`, .) %>%
    sweep_sparse(margin = 1, stats = ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat_list[[1]])) %>%
    magrittr::set_rownames(unique(ct_groups_n$ct))
  
  print("line 178")
  # indicator matrix for out group
  if (is.null(out_group_mat)) {
    out_group_mat <- matrix(1, nrow = dim(ct_groups_n)[1], ncol = dim(ct_groups_n)[1]) -
      diag(nrow = dim(ct_groups_n)[1], ncol = dim(ct_groups_n)[1])
  }else{
    #check if the outgroup is a binary matrix with all cell types are included in the rows and columns
    if (!is.matrix(out_group_mat) || any(out_group_mat != 0 & out_group_mat != 1)){
      stop("The out_group_mat must be a binary matrix with all cell types included in the rows and columns.")
    }
    if (! all(ct_groups_n$ct %in% rownames(out_group_mat)) || ! all(ct_groups_n$ct %in% colnames(out_group_mat))){
      stop("The out_group_mat must have all cell types included in the rows and columns.")
    }else{
      out_group_mat <- out_group_mat[ct_groups_n$ct, ct_groups_n$ct]
    }
  }
  
  #return(list(var_mat = var_mat, ct_groups_n = ct_groups_n, mean_mat = mean_mat, out_group_mat = out_group_mat, sum_mat = sum_mat))
  # mean for out group per cell type
  out_mean <- (out_group_mat %*% sum_mat) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% ct_groups_n$N, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat_list[[1]])) %>%
    magrittr::set_rownames(unique(ct_groups_n$ct))
  
  print("line 201")
  # variance for out groups per cell type
  #return(list(var_mat = var_mat, ct_groups_n = ct_groups_n, mean_mat = mean_mat, out_group_mat = out_group_mat, sum_mat = sum_mat, out_mean = out_mean))
  
  tot_var <- sweep_sparse(var_mat, margin = 1, stats = ct_groups_n$N - 1, fun = "*")
  print("line 207")
  tot_mean_sq <- sweep_sparse(mean_mat^2, margin = 1, stats = ct_groups_n$N, fun = "*")
  out_variance <- (out_group_mat %*% tot_var + out_group_mat %*% tot_mean_sq -
    sweep_sparse(out_mean^2, margin = 1, stats = out_group_mat %*% ct_groups_n$N, fun = "*")) %>%
    sweep_sparse(margin = 1, stats = out_group_mat %*% ct_groups_n$N - 1, fun = "/") %>%
    magrittr::set_colnames(rownames(data_mat_list[[1]])) %>%
    magrittr::set_rownames(unique(ct_groups_n$ct))
  print("line 210")
  # probability of relatively higher expression for each gene in each cell type (vs other genes)
  rel_exp <- (mean_mat - out_mean) /
    sqrt(sweep_sparse(var_mat, margin = 1, stats = ct_groups_n$N - 1, fun = "/") +
      sweep_sparse(out_variance, margin = 1, stats = out_group_mat %*% ct_groups_n$N - 1, fun = "/"))
  print("line 215")
  # calculate seismic specificity score
  #return(list(rel_exp = rel_exp, ratio_mat = ratio_mat))
  spec_score <- stats::pnorm(as.matrix(rel_exp)) * as.matrix(ratio_mat)
  
  out_spec_score <- out_group_mat %*% spec_score
  spec_score <- spec_score / (spec_score + out_spec_score)
  #spec_score <- sweep_sparse(x = spec_score, margin = 2, stats = Matrix::colSums(spec_score), fun = "/")

  return(list(Matrix::t(spec_score), Matrix::t(mean_mat), stats.dt))
}
