#' Compute group level statistics for SingleCellExperiment
#'
#' @param data_obj SingleCellExperiment object
#' @param meta_data Metadata dataframe. Rownames of the data frame should be the same as the colnames of the expression matrix (keep cell ID order correct).
#' @param group The column name in metadata for grouping
#' @param assay_name The name of the assay to extract
#' @param output_type The output type: either 'sce' or 'tibble'
#' @param filter_rare_ct If cell groups with too few cells should be filtered out. By default this is set to TRUE because these groups may be contaminants or something else. 
#' @param filter_thres If filter_rare_ct is TRUE, what cells we are going to filter out?
#' @param mean_only If only mean expression value calculated?
#' @return A SingleCellExperiment object or a tibble, depending on output_type
#' @export
#' 
cal_stat = function(data_obj,  group, meta_data = NULL, assay_name = "logcounts", output_type = "sce", filter_rare_ct = TRUE, filter_thres = 20, mean_only=FALSE) {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  
  #check assay_name exists in data_obj
  if (!assay_name %in% SummarizedExperiment::assayNames(data_obj)) {
    stop(paste0("Assay '", assay_name, "' does not exist in the SingleCellExperiment object."))
  }
  
  #check dimension match
  if (is.null(meta_data)){
    if(!group %in% colnames(SummarizedExperiment::colData(data_obj))){
      stop("Dimensions don't match or 'idents' variable is wrong. Check your input.")
    }
  }else if( (!group %in% colnames(meta_data) ) | (ncol(data_obj) != nrow(meta_data)) ){
    stop("Dimensions don't match or 'idents' variable is wrong. Check your input.")
  }
  
  #check name match for cells 
  if (!is.null(colnames(data_obj)) & !is.null(meta_data) & any(colnames(data_obj) != rownames(meta_data))){
    message("There seems to be a mismatch between cells in the data_mat and the meta_data. Check it carefully.")
  }
  
  data_mat = SummarizedExperiment::assay(data_obj, assay_name)
  feature_name = rownames(SummarizedExperiment::rowData(data_obj))
  #calculate cell num
  if (is.null(meta_data)){
    ident_group = SummarizedExperiment::colData(data_obj)[[group]]
    if(any(is.na(ident_group))){
      ident_group[which(is.na(ident_group))] = "NA"
      message("The group metadata seem to include NA ones, which have been replaced with strings instead.")
    }
  }else{
    ident_group = meta_data %>% dplyr::pull(group)
  }
  factor_mat = Matrix::fac2sparse(factor(ident_group, levels = unique(ident_group))) #model matrix
  cell_num = factor_mat %*% rep(1, length(ident_group)) %>% as.vector() %>% stats::setNames(unique(ident_group))
  #calculate mean
  #sum_mat = factor_mat %*% Matrix::t(data_mat)
  sum_mat = Matrix::t(data_mat %*% Matrix::t(factor_mat) )
  mean_mat = sum_mat %>% sweep_sparse(margin = 1, stats = cell_num , fun="/") %>%  #get mean expression matrix
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  
  ##filter cell groups
  if(!filter_rare_ct){
    filter_thres = 2 #only drop the singleton
  }
  ct_to_keep_idx = which(cell_num>= filter_thres)
  
  if(mean_only){ #if variance and expression ratio are not to computed
    cell_num = cell_num[ct_to_keep_idx]
    mean_mat = mean_mat[ct_to_keep_idx,]
    if (output_type=="tibble"){
      out_stat = cell_num %>%
        dplyr::as_tibble(rownames = "cellgroup") %>%
        dplyr::rename(cellnum = value) %>%
        dplyr::left_join(mean_mat %>% as.matrix %>% dplyr::as_tibble(rownames="cellgroup") %>% tidyr::pivot_longer(cols=all_of(feature_name), names_to = "genename", values_to="mean"), by=c("cellgroup"="cellgroup"), multiple = "all")
      return(out_stat)
    }else{
      #S4Vectors::metadata(data_obj)[["group_info"]] = list(cell_num, mean_mat) %>% purrr::set_names(c("cell_num","mean_mat"))
      data_obj = list(cell_num, mean_mat) %>% 
        purrr::set_names(c("cell_num","mean_mat")) %>%
        set_meta_slot(data_obj,"group_info",value=.)
    }
    return(data_obj)
  }
  
  #calculate variance 
  #var_mat = (factor_mat %*% (data_mat - t(factor_mat) %*% mean_mat)^2) %>%    #get variance matrix
  #var_mat =  factor_mat %*% Matrix::t(data_mat^2) 
  var_mat = (Matrix::t(data_mat^2 %*% Matrix::t(factor_mat) ) - 2*mean_mat*sum_mat + sweep_sparse(mean_mat^2, margin = 1, stats = cell_num , fun="*")) %>%  
    sweep_sparse(margin = 1, stats = cell_num-1 , fun="/") %>% 
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  #calculate ratio 
  #ratio_mat = (factor_mat %*% (Matrix::t(data_mat) > 0))%>% 
  ratio_mat = Matrix::t( (data_mat >0) %*% Matrix::t(factor_mat) ) %>% 
    sweep_sparse(margin = 1, stats = cell_num , fun="/") %>%
    magrittr::set_colnames(feature_name) %>% 
    magrittr::set_rownames(unique(ident_group ))
  
  var_mat = var_mat[ct_to_keep_idx,]
  ratio_mat = ratio_mat[ct_to_keep_idx,]
  cell_num = cell_num[ct_to_keep_idx]
  mean_mat = mean_mat[ct_to_keep_idx,]
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
    #set meta info
    if(meta_slot_is_null(data_obj, "obj_log")){
      obj_log_list = vector("list",length=7) %>%
        magrittr::set_names(c("group","progress","assay_name","rare_cell_filter_thres","out_group_mat","gene_trans_info","asso_model"))
      obj_log_list[["asso_model"]] =  vector("list",length=2) %>%
        magrittr::set_names(c("linear","spearman"))
    }else{
      obj_log_list = get_meta_slot(data_obj,"obj_log")
      obj_log_list[["gene_trans_info"]] = NULL
    }
    obj_log_list[["group"]] = group
    obj_log_list[["progress"]] = "cal_stat()"
    obj_log_list[["assay_name"]] = assay_name
    obj_log_list[["rare_cell_filter_thres"]] = filter_thres
  
    data_obj = list(cell_num, mean_mat, var_mat, ratio_mat) %>%
      purrr::set_names(c("cell_num","mean_mat","var_mat","ratio_mat")) %>%
      set_meta_slot(data_obj,"group_info",value=.) %>%
      set_meta_slot("obj_log",value = obj_log_list)
  }
  return(data_obj)
}
