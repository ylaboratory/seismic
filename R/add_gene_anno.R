#' Add gene-level statistics or annotation to the SingleCellExperiment object
#' The stored elements is in metadata (seismicGWAS.data) as a data frame (named 'gene_anno'). 
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
  #if (!meta_slot_is_null(data_obj,"gene_info")){
  if (!seismic_slot_is_null(data_obj,"gene_info")){
    #common_feature = intersect(feature_name, S4Vectors::metadata(data_obj)[["gene_info"]][["gene_name"]]) 
    common_feature = intersect(feature_name, get_seismic_slot(data_obj,"gene_info")[["gene_name"]]) 
    if (length(common_feature)==0){
      stop("Something seemed to be wrong, your gene index does not match the existing one")
    }else if(length(common_feature) < 0.5* length(feature_name)){
      warning("The gene ids of the gene annotation data (gene_info in seismicGWAS.data) poorly map the existing one.")
    }
  }
  
  ## add annotation
  #if (meta_slot_is_null(data_obj,"gene_info")){
  if (seismic_slot_is_null(data_obj,"gene_info")){
    #data_obj = set_meta_slot(data_obj, "gene_info", data.frame(gene_name = feature_name))
    data_obj = set_seismic_slot(data_obj, "gene_info", data.frame(gene_name = feature_name))
  }
  if (is.vector(gene_anno)){
    name = deparse(substitute(gene_anno))
    gene_anno = data.frame(gene_name = feature_name, value = gene_anno) %>% dplyr::rename_with(~name, .cols=2)
  }
  
  #remove duplicate columns in gene_info: update them
  #duplicated_col = colnames(get_meta_slot(data_obj,"gene_info"))[-1] %>% intersect(colnames(gene_anno))  
  duplicated_col = colnames(get_seismic_slot(data_obj,"gene_info"))[-1] %>% intersect(colnames(gene_anno))  
  if (length(duplicated_col) >=1 ){
    #data_obj = set_meta_slot(data_obj,"gene_info", get_meta_slot(data_obj,"gene_info") %>% dplyr::select(-any_of(duplicated_col)))
    data_obj = set_seismic_slot(data_obj,"gene_info", get_seismic_slot(data_obj,"gene_info") %>% dplyr::select(-any_of(duplicated_col)))
  }
  
  if (is.null(match_col)){
    #data_obj = get_meta_slot(data_obj,"gene_info") %>%
    data_obj = get_seismic_slot(data_obj,"gene_info") %>%
      dplyr::full_join(gene_anno, by="gene_name") %>%
      as.data.frame() %>%
      #set_meta_slot(data_obj,slot="gene_info",value=.)
      set_seismic_slot(data_obj,slot="gene_info",value=.)
  }else{
    #data_obj = get_meta_slot(data_obj,"gene_info") %>%
    data_obj = get_seismic_slot(data_obj,"gene_info") %>%
      dplyr::full_join(gene_anno, by=c("gene_name"=match_col)) %>%
      as.data.frame() %>%
      #set_meta_slot(data_obj,slot="gene_info",value=.)
      set_seismic_slot(data_obj,slot="gene_info",value=.)
  }
  
  
  return(data_obj)
}

