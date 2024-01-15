#' Add some specific global statistics to the metadata 'gene_anno' slot. 
#' 
#' Some predefined statistics are:
#' 'det_cell_num': number of cells that express such a gene
#' 'ave_exp_ct': mean value of a gene's average expression across all cell types 
#' 'max_exp_ct': maximum value of a gene's average expression across all cell types 
#' 'ave_rat_ct': mean value of the ratio of cells that express the gene across all cell types
#' 'max_rat_ct': maximum value of the ratio of cells that express the gene across all cell types
#' 'ave_exp_all': mean expression value of a gene across all cells (the data set)
#' 
#' @param data_obj SingleCellExperiment object
#' @param stats A vector that select the predefined statistics for the current data set.
#' @return A SingleCellExperiment object or a tibble. With meta data slots are translate into the 
#' @export
#' 
add_glob_stats = function(data_obj, stats=NULL) {
  if (meta_slot_is_null(data_obj,"group_info")){
    stop("Please run cal_stat() first")
  }
  if (is.null(stats) | (!all(stats %in% c("det_cell_num","ave_exp_ct","max_exp_ct","ave_rat_ct","max_rat_ct","ave_exp_all")))){
    stop("Please correctly specify what statistics you are going to calculate")
  }
  
  if ("det_cell_num" %in% stats){
    det_cell_num = sweep_sparse(get_meta_slot(data_obj,"group_info")[["ratio_mat"]],margin=1, stats = get_meta_slot(data_obj,"group_info")[["cell_num"]], fun = "*")
    det_cell_num = Matrix::colSums(det_cell_num)
    data_obj = add_gene_anno(data_obj, gene_anno = det_cell_num)
  }
  
  if ("ave_exp_ct" %in% stats){
    ave_exp_ct = Matrix::colMeans(get_meta_slot(data_obj,"group_info")[["mean_mat"]])
    data_obj = add_gene_anno(data_obj, gene_anno = ave_exp_ct)
  }
  
  if ("max_exp_ct" %in% stats){
    max_exp_ct = sparseMatrixStats::colMaxs(get_meta_slot(data_obj,"group_info")[["mean_mat"]]) %>% purrr::set_names(colnames(get_meta_slot(data_obj,"group_info")[["mean_mat"]]))
    data_obj = add_gene_anno(data_obj, gene_anno = max_exp_ct)
  }
  
  if ("ave_rat_ct" %in% stats){
    ave_rat_ct = Matrix::colMeans(get_meta_slot(data_obj,"group_info")[["ratio_mat"]])
    data_obj = add_gene_anno(data_obj, gene_anno = ave_rat_ct)
  }
  
  if ("max_rat_ct" %in% stats){
    max_rat_ct = sparseMatrixStats::colMaxs(get_meta_slot(data_obj,"group_info")[["ratio_mat"]])
    data_obj = add_gene_anno(data_obj, gene_anno = max_rat_ct)
  }
  
  if ("ave_exp_all" %in% stats){
    ave_exp_all = sweep_sparse(get_meta_slot(data_obj,"group_info")[["mean_mat"]],margin=1, stats = get_meta_slot(data_obj,"group_info")[["cell_num"]]) %>%
      Matrix::colSums() %>%
      `/`(sum(get_meta_slot(data_obj,"group_info")[["cell_num"]])) 
    data_obj = add_gene_anno(data_obj, gene_anno = ave_exp_all)
  }
  
  return(data_obj)
}