##seismicGWAS interface and functions for better access the metadata slot of SCE

#' Check if a metadata slot is empty or not.
#' @param data_obj SingleCellExperiment object
#' @param slot slot to retrieve in metadata (where seismicGWAS stores the information), any one of "group_info","gene_info","association" or ""
#' @return TRUE or FALSE
#' @keywords internal
#' @noRd
meta_slot_is_null = function(data_obj, slot){
  if (!slot %in% c("group_info","gene_info","association","cell_type_anno","obj_log")){
    stop("slot is wrong")
  }
  return(is.null(S4Vectors::metadata(data_obj)[[slot]]))
}


#' Get the value of a metadata slot
#' @param data_obj SingleCellExperiment object
#' @param slot  slot to retrieve in metadata (where seismicGWAS stores the information), any one of "group_info","gene_info","association"
#' @return the value of the sloot
#' @keywords internal
#' @noRd
get_meta_slot = function(data_obj, slot){
  if (!slot %in% c("group_info","gene_info","association","cell_type_anno","obj_log") ){
    stop("slot is wrong")
  }
  return(S4Vectors::metadata(data_obj)[[slot]])
}

#' Set the value of a metadata slot
#' @param data_obj SingleCellExperiment object
#' @param slot  slot to retrieve in metadata (where seismicGWAS stores the information), any one of "group_info","gene_info","association"
#' @param value Value of the new slot
#' @return the data object with a new set value
#' @keywords internal
#' @noRd
set_meta_slot = function(data_obj, slot, value){
  if (!slot %in% c("group_info","gene_info","association","cell_type_anno","obj_log")){
    stop("slot is wrong")
  }
  S4Vectors::metadata(data_obj)[[slot]] = value
  return(data_obj)
}


#' add a cell type association data frame to data_obj
#' @param data_obj SingleCellExperiment object
#' @param ct_asso_df A list or a single data frame. If it's a single data frame then requires a name else it requires names for each list. 
#' The data frame has to have two columns: cell_type and Pvalue 
#' @param trait_name name of the trait(s)
#' @param asso_model The model used to calculate the cell type level association with the trait.
#' @return the data object with set values
#' @keywords internal
#' @noRd
add_ct_asso = function(data_obj, ct_asso_df, trait_name, asso_model){ 
  if (!asso_model %in% c("linear","spearman")){
    stop("Please indicate the right model.")
  }
  full_asso_list =  get_meta_slot(data_obj,"association")

  all_asso_df = full_asso_list[[asso_model]] %>% 
    keep(~!is.null(names(.)))
  all_trait_name = c(names(all_asso_df), trait_name)
  if(is.null(all_asso_df) & !is.list(ct_asso_df)){
    all_asso_df = list(ct_asso_df)
  }else{
    all_asso_df = all_asso_df %>% append(ct_asso_df) 
  }
  all_asso_df = all_asso_df %>%
    purrr::set_names(all_trait_name) %>% #set names
    .[!duplicated(names(.), fromLast=TRUE)] #remove duplicated traits
  
  full_asso_list[[asso_model]] = all_asso_df 
  full_asso_list = full_asso_list[c("linear","spearman")[names(full_asso_list) %in% c("linear","spearman")]] #rearrange 
  return(set_meta_slot(data_obj,"association",full_asso_list))
}


#' Get the data frame (s) of cell type-level association created by seismicGWAS.
#' @param data_obj SingleCellExperiment object
#' @param trait_name Name of the trait to get the association values. It could also be a list of traits or "all".
#' @param asso_model The model previously specified too calculate the cell type level association with the trait.
#' @param merge_output Default is FALSE. But if you wish to get a single full data frame (instead of a list) of all P_value of a list of traits, this can be set into TRUE.
#' @return the data object with set values
#' @export 
get_ct_asso = function(data_obj, trait_name, asso_model, merge_output = FALSE){
  if(is.null(data_obj)|is.null(trait_name)|is.null(asso_model)){
    stop("Input should not be NULL")
  }
  if(!inherits(data_obj,"SingleCellExperiment")){
    stop("object class fault")
  }
  trait_name = unlist(trait_name)
  
  if(!asso_model %in% names(get_meta_slot(data_obj,"association")) ){
    stop("You have not calculated the cell type association using this model")
  }
  
  if(length(trait_name)==1){
    if(trait_name!="all"){
      if( !trait_name %in% names(get_meta_slot(data_obj,"association")[[asso_model]])){
        stop("The trait are not present in the cell type association tables")
      }else{
        return(get_meta_slot(data_obj,"association")[[asso_model]][[trait_name]])
      }
    }else{
      asso_df_list = get_meta_slot(data_obj,"association")[[asso_model]]
    }
  }else{
    if (any(!trait_name %in% names(get_meta_slot(data_obj,"association")[[asso_model]]))){
      stop("Some traits are not present in the cell type association tables")
    }
    asso_df_list = get_meta_slot(data_obj,"association")[[asso_model]][trait_name]
  }
  #merge output or not 
  if(merge_output){
    asso_df_list = asso_df_list %>%
      purrr::map(~select(.x, cell_type,Pvalue)) %>%
      purrr::map2(names(.), ~magrittr::set_colnames(.x,c("cell_type",.y))) 
    #check if the group contradict
    all_cell_names = asso_df_list %>%
      purrr::map(~pull(.x,"cell_type")) 
    if(length(all_cell_names[[1]])!= length(unique(purrr::reduce(all_cell_names, ~c(.x,.y))))){
      warning("It seems like some cell type level association for some traits did not match. Probably the results are not from the same time of analysis. Merge output is not recommended")
    }
    #print out results
    asso_df_list = asso_df_list %>%
      purrr::reduce(~full_join(.x,.y,by="cell_type"))
  }
  return(asso_df_list)
}

####user-friendly interface to retrieve the information
#' Print out the summary information so far for the processing of seismicGWAS
#' @param data_obj SingleCellExperiment object
#' @param verbose How detailed the information is being printed? If this is set to be False, only part of the running log will be return.
#' The more detailed information will be printed out if this is set to True.
#' @param info_to_return If true, the whole metadata
#' @return The returned value depends on the input. If info_to_print parameter is null, there will be no values or data returned.
#' @export 
seismic_summary_info = function(data_obj, verbose = F, info_to_return = F){
  if(meta_slot_is_null(data_obj,"obj_log")){
    message("You have not started analysis, no log information exist")
    return()
  }
  obj_log_list = get_meta_slot(data_obj,"obj_log")
  message("The data set includes information of ",nrow(facs.sce)," genes and ",ncol(facs.sce),"cells. \n")
  if(verbose){
    message("The genes include ",paste(head(rownames(facs.sce)),collapse = ", "), ".\n")
  }
  message("And the analysis granularity is",obj_log_list[["group"]], " in ",obj_log_list[["assay_name"]]," assay \n")
  if(verbose){
    if( !is.null(obj_log_list[["rare_cell_filter_thres"]])){
      message("After filtering of rare cell types with fewer than ",obj_log_list[["rare_cell_filter_thres"]]," cells, there are still ",sum(get_meta_slot(data_obj,"group_info")[["cell_num"]])," cells left. \n")
    }
    message("The cell types include ",paste(head(names(get_meta_slot(data_obj,"group_info")[["cell_num"]])),collapse = ", "), ".\n")
  }
  message(paste0("The previously step is ",obj_log_list[["progress"]]," \n"))
  if(obj_log_list[["progress"]]!="cal_stat()"){
    if(is.null(obj_log_list[["out_group_mat"]])){
      message("Basic specificity score is calculated. \n")
    }else{
      message("A out group matrix was used to specify the out group assignment for different cell types. \n")
    }
    if(verbose){
      message("Genes in the specificity score matrix now include ",paste(head(colnames(t)),collapse = ", "),".\n")
    }
  }
  #if has calculated the association
  if(obj_log_list[["progress"]]=="cal_ct_asso()" ){
    message("Cell type associations have been calculated using the ",paste(names(obj_log_list[["asso_model"]] %>% purrr::keep(~is.null(.x))),collapse = ", ")," model(s). \n")
    for (model in obj_log_list[["aaso_model"]]){
      if(!is.null(model)){
        message("Association have been calculated using ",length(model[["model_genes"]])," genes with the ",model," model in ", length(model[["traits"]])," traits.")
        if(verbose){
          message("The traits include: ", paste(model[["traits"]], collapse = ", "),".")
        }
        message("\n")
      }
    }
  } 
  
  if(info_to_return){
    return(obj_log_list)
  }
}