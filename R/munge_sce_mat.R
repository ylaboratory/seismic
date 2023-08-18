#' Map feature (gene) names (row names) to either entrez or other format
#' 
#' This function is similar to `trans_mmu_to_hsa_stat`.
#' A new SingleCellExperiment object will be created, with cell meta data the same as the original 
#' The gene meta data will be discarded
#'
#' @param data_obj SingleCellExperiment object
#' @param assay_name name of the assay for the mapping. By default it will be all assays. 
#' @param mapping_df A data frame to specify the mapping. The first column should contain all the feature names while the second column indicates the new feature names to be mapped to.
#' @param multi_mapping How to translate feature name when multiple mapping exists?? By "mean" value or "sum" value? (sum for default)
#' @return A SingleCellExperiment object 
#' @importFrom S4Vectors DataFrame
#' @export
#' 
munge_sce_mat = function(data_obj,  mapping_df, assay_name = "all",multi_mapping = "sum") {
  #check if the assay exists
  if ( assay_name != "all" & !assay_name %in% SummarizedExperiment::assayNames(data_obj)){
    stop("The assay you are indicating does not exist")
  }
  #check if the feature name is correct
  if( is.null(rownames(data_obj)) | !any( rownames(data_obj) %in% (mapping_df %>% dplyr::pull(1))) ){
    stop("The feature names do not match the first column of the mapping_df")
  }
  if (assay_name=="all"){
    assay_name =  SummarizedExperiment::assayNames(data_obj)
  }
  
  #merge
  new_assay = list()
  mapping_df = mapping_df %>% mutate_all(~as.character(.))
  
  for (assay_name_i in assay_name){
    data_mat = SummarizedExperiment::assay(data_obj, assay_name_i)
    #filter mapping
    all_mapping = mapping_df %>% 
      dplyr::filter(.[[1]] %in% rownames(data_mat)) %>% #filter by feature names
      tidyr::drop_na(2) #drop features without mapping
    #subset rows
    data_mat = data_mat[match(all_mapping[[1]],rownames(data_mat)),]
    #split matrix into two
    multiple_mapping = all_mapping %>% 
      group_by_at(2) %>%
      filter(n()>=2) %>%
      ungroup
    #separate the matrix with multiple mapping a
    single_mat = data_mat[which(!rownames(data_mat) %in% multiple_mapping[[1]]),]
    multi_mat = data_mat[which(rownames(data_mat) %in% multiple_mapping[[1]]),]
    
    #transform mapping 
    fac_mat = Matrix::fac2sparse(factor(multiple_mapping[[2]], levels = unique(multiple_mapping[[2]]))) 
    if(multi_mapping=="mean"){
      fac_mat = sweep_sparse(fac_mat, margin=1, stats = Matrix::rowSums(fac_mat), fun = "/")
    }
    #merge and set rownames
    multi_mat = ( fac_mat %*% multi_mat) %>% 
      magrittr::set_rownames(unique(multiple_mapping[[2]])) %>%
      magrittr::set_colnames(colnames(data_mat))
    single_mat = single_mat %>% magrittr::set_rownames(all_mapping %>% filter(! .[[1]] %in% multiple_mapping[[1]] ) %>% pull(2))
    data_mat = rbind(single_mat, multi_mat)
    #rearrange based on the new gene names
    data_mat = data_mat[match(unique(all_mapping[[2]]),rownames(data_mat)),]
    #add assay
    new_assay = new_assay %>% append(data_mat)
  }
  new_assay = new_assay %>% purrr::set_names(assay_name)
  
  #create output data object
  out_data = SingleCellExperiment::SingleCellExperiment(
    assays = new_assay,
    colData = SummarizedExperiment::colData(data_obj),
    rowData = data.frame(rownames(new_assay[[1]]), row.names = rownames(new_assay[[1]]),fix.empty.names = FALSE) %>% 
      magrittr::set_colnames(colnames(mapping_df)[2]) %>%
      S4Vectors::DataFrame()
  )
  
  return(out_data)
}


