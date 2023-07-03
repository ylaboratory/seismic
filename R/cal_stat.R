#' Compute group level statistics for SingleCellExperiment
#'
#' @param data_obj SingleCellExperiment object
#' @param meta_data Metadata dataframe
#' @param group The column name in metadata for grouping
#' @param assay_name The name of the assay to extract
#' @param output_type The output type: either 'sce' or 'tibble'
#' @param filter_rare_ct If cell groups with too few cells should be filtered out. By default this is set to TRUE because these groups may be contaminants or something else. 
#' @param filter_thres If filter_rare_ct is TRUE, what cells we are going to filter out?
#' @return A SingleCellExperiment object or a tibble, depending on output_type
#' @export
#' 
cal_stat = function(data_obj, meta_data = NULL, group, assay_name = "logcounts", output_type = "sce", filter_rare_ct = TRUE, filter_thres = 20) {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  
  #check assay_name exists in data_obj
  if (!assay_name %in% assayNames(data_obj)) {
    stop(paste0("Assay '", assay_name, "' does not exist in the SingleCellExperiment object."))
  }
  
  #check dimension match
  if((is.null(meta_data) & !group %in% colnames(colData(data_obj)) )| (! group %in% colnames(meta_data) )| (ncol(data_obj) != nrow(meta_data))){
      stop("Dimensions don't match or 'idents' variable is wrong. Check your input.")
  }
  
  #check name match for cells 
  if (!is.null(colnames(data_obj)) & !is.null(meta_data) & any(colnames(data_obj) != rownames(meta_data))){
    message("There seems to be a mismatch between cells in the data_mat and the meta_data. Check it carefully.")
  }
  
  data_mat = assay(data_obj, assay_name)
  feature_name = rownames(rowData(data_obj))
  #calculate cell num
  if (is.null(meta_data)){
    ident_group = colData(data_obj)[[group]]
  }else{
    ident_group = meta_data %>% dplyr::pull(group)
  }
  factor_mat = Matrix::fac2sparse(factor(ident_group, levels = unique(ident_group))) #model matrix
  cell_num = factor_mat %*% rep(1, length(ident_group)) %>% as.vector() %>% stats::setNames(unique(ident_group))
  #calculate mean
  sum_mat = factor_mat %*% Matrix::t(data_mat)
  mean_mat = sum_mat %>% sweep_sparse(margin = 1, stats = cell_num , fun="/") %>%  #get mean expression matrix
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  #calculate variance 
  #var_mat = (factor_mat %*% (data_mat - t(factor_mat) %*% mean_mat)^2) %>%    #get variance matrix
  var_mat = (factor_mat %*% Matrix::t(data_mat)^2 - 2*mean_mat*sum_mat + sweep_sparse(mean_mat^2, margin = 1, stats = cell_num , fun="*")) %>%  
    sweep_sparse(margin = 1, stats = cell_num-1 , fun="/") %>% 
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  #calculate ratio  
  ratio_mat = (factor_mat %*% (Matrix::t(data_mat) > 0))%>% 
    sweep_sparse(margin = 1, stats = cell_num , fun="/") %>%
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  
  ##filter cell groups
  if(!filter_rare_ct){
    filter_thres = 1 #only drop the singleton
  }
  ct_to_keep_idx = which(cell_num>= filter_thres)
  cell_num = cell_num[ct_to_keep_idx]
  mean_mat = mean_mat[ct_to_keep_idx,]
  var_mat = var_mat[ct_to_keep_idx,]
  ratio_mat = ratio_mat[ct_to_keep_idx,]
  
  ##output  
  if (output_type=="tibble"){
    out_stat = cell_num %>%
      dplyr::as_tibble(rownames = "cellgroup") %>%
      dplyr::rename(cellnum = value) %>%
      dplyr::left_join(mean_mat %>% as.matrix %>% dplyr::as_tibble(rownames="cellgroup") %>% tidyr::pivot_longer(cols=all_of(feature_name), names_to = "genename", values_to="mean"), by=c("cellgroup"="cellgroup"), multiple = "all") %>%
      dplyr::left_join(var_mat %>% as.matrix %>% dplyr::as_tibble(rownames="cellgroup") %>% tidyr::pivot_longer(cols=all_of(feature_name), names_to = "genename", values_to="variance"), by=c("cellgroup"="cellgroup","genename"="genename")) %>%
      dplyr::left_join(ratio_mat %>% as.matrix %>% dplyr::as_tibble(rownames="cellgroup") %>% tidyr::pivot_longer(cols=all_of(feature_name), names_to = "genename", values_to="ratio"), by=c("cellgroup"="cellgroup","genename"="genename")) 
    return(out_stat)
  }else{
    metadata(data_obj)[["group_info"]] = list(cell_num, mean_mat, var_mat, ratio_mat) %>% purrr::set_names(c("cell_num","mean_mat","var_mat","ratio_mat"))
  }
  return(data_obj)
}
