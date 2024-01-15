#' Compute group level specificity score for genes.
#' This can be enhanced by specifying the rule to choose out groups (cell groups that are considered as background cells for a cell group)
#'
#' @param data_obj SingleCellExperiment object or a tibble from cal_stat() function. The meta data contain the basic statistics for specificity score computation.
#' @param out_group_mat A matrix to specify how to choose out groups. Each row indicate if all other cell groups (in the corresponding column) are (TRUE) or are not (FALSE) out groups. 
#'                      Rows and columns should be named with the cell group names in the previous step. It's recommended to add both row names and column names to the out group matrix.
#'                      If there aren't any names, the rows and columns should present the same order as in the cell groups (given in S4Vectors::metadata(data_obj)[["group_info"]])
#' @return A SingleCellExperiment object or a tibble, depending on the input type
#' @export

cal_sscore = function(data_obj, out_group_mat = NULL){
  #check if the meta data is correct
  #if (is.null(S4Vectors::metadata(data_obj)[["group_info"]]) | any(!c("cell_num","mean_mat","var_mat","ratio_mat") %in% names(S4Vectors::metadata(data_obj)[["group_info"]]))){
  if (meta_slot_is_null(data_obj,"group_info") | any(!c("cell_num","mean_mat","var_mat","ratio_mat") %in% names(get_meta_slot(data_obj,"group_info")))){
    stop("Please run cal_stat() first")
  }
  
  #reshuffle out group matrix and check if it's correct
  #ct_names = S4Vectors::metadata(data_obj)[["group_info"]][["cell_num"]] %>% names
  ct_names = names(get_meta_slot(data_obj,"group_info")[["cell_num"]])
  if(!is.null(out_group_mat)){
    if((!is.null(rownames(out_group_mat)) & any(! ct_names %in% rownames(out_group_mat))) | (is.null(rownames(out_group_mat)) & nrow(out_group_mat)!=length(ct_names) )){
      stop("The setting of rownames is wrong.")
    }else if (!is.null(rownames(out_group_mat))){
      out_group_mat = out_group_mat[match(ct_names, rownames(out_group_mat)),]
    }
    if ((!is.null(colnames(out_group_mat)) & any(! ct_names %in% colnames(out_group_mat))) | (is.null(colnames(out_group_mat)) & ncol(out_group_mat)!=length(ct_names) )){
      stop("The setting of column names is wrong.")
    }else if (!is.null(colnames(out_group_mat))){
      out_group_mat = out_group_mat[,match(ct_names, colnames(out_group_mat))]
    }
  }else{
    out_group_mat = matrix(1,nrow = length(ct_names), ncol = length(ct_names)) - diag(nrow = length(ct_names), ncol = length(ct_names))
  }
  
  #calculate relative expression
  cell_num = get_meta_slot(data_obj,"group_info")[["cell_num"]]
  mean_mat = get_meta_slot(data_obj,"group_info")[["mean_mat"]]
  var_mat = get_meta_slot(data_obj,"group_info")[["var_mat"]]
  ratio_mat = get_meta_slot(data_obj,"group_info")[["ratio_mat"]]
  
  ##out_group_num = sum(cell_num) - cell_num
  #all_ones = matrix(1, nrow = length(cell_num), ncol = length(cell_num))
  
  #calculate out mean 
  #tot_mean = sweep_sparse(mean_mat, margin = 1, stats = cell_num, fun="*")
  #out_mean =  (all_ones %*% tot_mean - tot_mean) %>% sweep_sparse(margin = 1, stats = out_group_num, fun="/")
  tot_mean = sweep_sparse(mean_mat, margin = 1, stats = cell_num, fun="*")
  out_mean = (out_group_mat %*% tot_mean) %>% sweep_sparse(margin = 1, stats = out_group_mat %*% cell_num, fun="/")
  
  ##calculate out variance 
  #tot_var = sweep_sparse(var_mat, margin = 1, stats = cell_num-1, fun="*")
  #tot_mean_sq = sweep_sparse(mean_mat^2, margin = 1, stats = cell_num, fun="*")
  #out_variance = (all_ones %*% tot_var - tot_var + all_ones %*% tot_mean_sq - tot_mean_sq - sweep_sparse(out_mean^2, margin=1, stat = out_group_num, fun="*")) %>% 
  #tot_var = sweep_sparse(margin=1,stats=out_group_num-1, fun="/")
  tot_var = sweep_sparse(var_mat, margin = 1, stats = cell_num-1, fun="*")
  tot_mean_sq = sweep_sparse(mean_mat^2, margin = 1, stats = cell_num, fun="*")
  out_variance = (out_group_mat %*% tot_var + out_group_mat %*% tot_mean_sq - sweep_sparse(out_mean^2, margin=1, stat = out_group_mat %*% cell_num,fun="*") ) %>%
    sweep_sparse(margin=1,stats=out_group_mat %*% cell_num -1, fun="/")
  ##relative expression
  #rel_exp = (mean_mat - out_mean)/sqrt(sweep_sparse(var_mat,margin=1, stats=cell_num-1,fun="/" ) + sweep_sparse(out_variance,margin = 1,stats = out_group_num-1, fun="/"))
  rel_exp = (mean_mat - out_mean)/sqrt(sweep_sparse(var_mat,margin=1, stats=cell_num-1,fun="/" ) + sweep_sparse(out_variance,margin = 1,stats =  out_group_mat %*% cell_num-1, fun="/"))
  
  ##specificity score
  sscore = pnorm(as.matrix(rel_exp))*ratio_mat 
  #sscore = sscore/(all_ones %*% sscore)
  sscore = sweep_sparse(x= sscore, margin = 2, stats=Matrix::colSums(sscore), fun="/")
  #S4Vectors::metadata(data_obj)[["group_info"]][["sscore"]] = sscore
  data_obj = get_meta_slot(data_obj,"group_info") %>%
    append(sscore) %>% 
    purrr::set_names(c(names(get_meta_slot(data_obj,"group_info")), "sscore")) %>%
    set_meta_slot(data_obj,"group_info", value=.)
    
  return(data_obj)
}
