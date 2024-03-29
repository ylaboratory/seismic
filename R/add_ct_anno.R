#' Add cell-type level annotation
#' The stored elements is in metadata (seismicGWAS.data) as a data frame (named 'ct_anno'). 
#'
#' @param data_obj SingleCellExperiment object
#' @param ct_anno A data frame, tibble or a vector to annotate features. 
#' @param match_col To specify, only when a tibble or a data frame input is given, which column should match that of the feature names? 
#' @return A SincleCellExperiment object
#' @export
#' 
add_ct_anno = function(data_obj, ct_anno, match_col=NULL) {
  #check input type
  if (! (is.data.frame(ct_anno)|is.vector(ct_anno)| dplyr::is.tbl(ct_anno))){
    stop("Wrong input 'ct_anno' type")
  }
  #check if 'match_col' exist in the df
  if (is.data.frame(ct_anno) | dplyr::is.tbl(ct_anno) ){ 
    if (is.null(match_col) & is.data.frame(ct_anno)){
      feture_name = rownames (ct_anno)
    }else{
      if (! match_col %in%colnames(ct_anno) ){
        stop("Please indicate a valid column in ct_anno")
      }
      ct_name = ct_anno %>% dplyr::select(all_of(match_col)) %>% pull
    }
  }else{
    if(is.null(names(ct_anno) )){
      stop("As vector input for the function, the names should be indicated")
    }
    ct_name = names(ct_anno)
  }
  #check if current index conflicts with the previous index
  
  #if (!any(ct_name %in% names(get_meta_slot(data_obj,"group_info")[["cell_num"]]) )){
  if (!any(ct_name %in% names(get_seismic_slot(data_obj,"group_info")[["cell_num"]]) )){
    stop("Something seems to be wrong, your cell types do not match the existing data. Or group-level information is not added to the current data set.")
  }
  
  ## add 
  #if (meta_slot_is_null(data_obj,"cell_type_anno")){
  if (seismic_slot_is_null(data_obj,"cell_type_anno")){
    #if (!meta_slot_is_null(data_obj,"group_info")){
    if (!seismic_slot_is_null(data_obj,"group_info")){
      #data_obj = set_meta_slot(data_obj,"cell_type_anno", data.frame(cell_type = names(get_meta_slot(data_obj,"group_info")[["cell_num"]])))
      data_obj = set_seismic_slot(data_obj,"cell_type_anno", data.frame(cell_type = names(get_seismic_slot(data_obj,"group_info")[["cell_num"]])))
    }else{
      #data_obj = set_meta_slot(data_obj,"cell_type_anno",data.frame(cell_type = ct_name))
      data_obj = set_seismic_slot(data_obj,"cell_type_anno",data.frame(cell_type = ct_name))
    }
  }
  if (is.vector(ct_anno)){
    name = rlang::as_name(rlang::ensym(ct_anno))
    ct_anno = data.frame(cell_type = ct_name, value = ct_anno) %>% dplyr::rename_with(~name, .cols=2)
  }
  
  #remove duplicate columns in gene_info: update them
  #duplicated_col = colnames(S4Vectors::metadata(data_obj)[["cell_type_anno"]])[-1] %>% intersect(colnames(ct_anno))  
  duplicated_col = colnames(get_seismic_slot(data_obj,"cell_type_anno"))[-1] %>% intersect(colnames(ct_anno))  
  if (length(duplicated_col) >=1 ){
    #data_obj = get_meta_slot(data_obj,"cell_type_anno") %>%
    data_obj = get_seismic_slot(data_obj,"cell_type_anno") %>%
      dplyr::select(-any_of(duplicated_col)) %>%
      #set_meta_slot(data_obj, "cell_type_anno",value=.)
      set_seismic_slot(data_obj, "cell_type_anno",value=.)
  }
  
  ##join
  if(is.null(match_col)){
    #data_obj = get_meta_slot(data_obj, "cell_type_anno") %>%
    data_obj = get_seismic_slot(data_obj, "cell_type_anno") %>%
      dplyr::full_join(ct_anno, by="cell_type") %>% 
      as.data.frame() %>%
      #set_meta_slot(data_obj, "cell_type_anno", value=.)
      set_seismic_slot(data_obj, "cell_type_anno", value=.)
  }else{
    #data_obj = get_meta_slot(data_obj,"cell_type_anno") %>%
    data_obj = get_seismic_slot(data_obj,"cell_type_anno") %>%
      dplyr::full_join(ct_anno, by=c("cell_type"=match_col)) %>%
      as.data.frame() %>%
      #set_meta_slot(data_obj, "cell_type_anno", value=.)
      set_seismic_slot(data_obj, "cell_type_anno", value=.)
  }
  
  return(data_obj)
}
