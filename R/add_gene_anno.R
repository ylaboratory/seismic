#' Add gene-level statistics or annotation to the SingleCellExperiment object
#' The stored elements is in meta data as a data frame (named 'gene_anno'). 
#'
#' @param data_obj SingleCellExperiment object
#' @param gene_anno A data frame, tibble or a vector to annotate features. 
#' @param match_col To specify, only when a tibble or a data frame input is given, which column should match that of the feature names? 
#' @return A SincleCellExperiment object,
#' @export
#' 
add_gene_anno = function(data_obj, gene_anno, match_col=NULL) {
  #check input type
  if (! (is.data.frame(gene_anno)|is.vector(gene_anno)| dplyr::is.tbl(gene_anno))){
    stop("Wrong input 'gene_anno' type")
  }
  #check if 'match_col' exist in the df
  if (is.data.frame(gene_anno) | dplyr::is.tbl(gene_anno) ){ 
    if (is.null(match_col) & is.data.frame(gene_anno)){
      feture_name = rownames (gene_anno)
    }else{
      if (! match_col %in%colnames(gene_anno) ){
        stop("Please indicate a valid column in gene_anno")
      }
      feature_name = gene_anno %>% dplyr::select(all_of(match_col)) %>% pull
    }
  }else{
    if(is.null(names(gene_anno) )){
      stop("As vector input for the function, the names should be indicated")
    }
    feature_name = names(gene_anno)
  }
  #check if current index conflicts with the previous index
  
  if ((!is.null(metadata(data_obj)[["gene_info"]])  ) & (!any(feature_name %in% metadata(data_obj)[["gene_info"]][["gene_name"]])) ){
    stop("Something seemed to be wrong, your gene index does not match the existing one")
  }
  
  ## add 
  if (is.null(metadata(data_obj)[["gene_info"]])){
    metadata(data_obj)[["gene_info"]] = data.frame(gene_name = feature_name)
  }
  if (is.vector(gene_anno)){
    name = deparse(substitute(gene_anno))
    gene_anno = data.frame(gene_name = feature_name, value = gene_anno) %>% dplyr::rename_with(~name, .cols=2)
  }
  
  #remove duplicate columns in gene_info: update them
  duplicated_col = colnames(metadata(data_obj)[["gene_info"]])[-1] %>% intersect(colnames(gene_anno))  
  if (length(duplicated_col) >=1 ){
    metadata(data_obj)[["gene_info"]] = metadata(data_obj)[["gene_info"]] %>% dplyr::select(-any_of(duplicated_col))
  }
  
  if (is.null(match_col)){
    metadata(data_obj)[["gene_info"]] = metadata(data_obj)[["gene_info"]] %>% dplyr::full_join(gene_anno, by="gene_name") %>% as.data.frame()
  }else{
    metadata(data_obj)[["gene_info"]] = metadata(data_obj)[["gene_info"]] %>% dplyr::full_join(gene_anno, by=c("gene_name"=match_col)) %>% as.data.frame()
  }
  
  return(data_obj)
}

