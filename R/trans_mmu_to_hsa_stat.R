#' Translate metrics named after mouse gene id to human gene id
#' Only genes with entrez ID are kept.
#' Mapping can be also done between symbol, ensembl id, etc
#' Features with no mapping will be discarded.
#'
#' @param data_obj SingleCellExperiment object
#' @param gene_mapping_table A data frame or tibble to represent gene mapping. 
#' @param from The column (as in gene_mapping_table) to translate the feature names from. 
#' @param to The column (as in gene_mapping_table) to translate the feature names to.
#' @param multi_mapping How to translate metrics when one feature is mapped to several? By "mean" value or "sum" value? (mean for default)
#' @param metadata_type Meta data slots you are translating. By default all meta data (all group-level statistics and specificity score) will be translate. You can also specify anyone if you only would like to translate it.
#' @return A SingleCellExperiment object or a tibble. With meta data slots are translate into the 
#' @export
#' 
trans_mmu_to_hsa_stat = function(data_obj, gene_mapping_table, from, to, multi_mapping = "mean", metadata_type = "all") {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  #check column names 
  if (!all(c(from,to) %in% colnames(gene_mapping_table))) {
    stop("The gene_mapping_table is not fine or the columns of 'from' and 'to' do not exist in the table")
  }
  #check metadata slots 
  if (is.null(get_meta_slot(data_obj,"group_info"))|(metadata_type!="all" & (!metadata_type %in% c("ratio","mean_mat","var_mat","ratio_mat","sscore")))){
    stop("The meta data slot does not exist")
  }
  #check feature names
  if (metadata_type=="all"){
    all_names = names(get_meta_slot(data_obj,"group_info")) %>%
      .[.!="cell_num"] %>%
      purrr::map(~colnames(get_meta_slot(data_obj,"group_info")[[.x]])) %>%
      purrr::reduce(~unique(c(.x,.y)))
  }else{
    all_names = colnames(get_meta_slot(data_obj,"group_info")[[metadata_type]])
  }
  if (!any(all_names %in% (gene_mapping_table %>% dplyr::pull(from)))){
    stop("Something's wrong. The feature names don't match the column in the gene_mapping_table")
  }
  #checek multi_mapping parameters
  if (!multi_mapping %in% c("mean","sum")) {
    stop("multi_mapping can only be either 'mean' or 'sum'")
  }
  
  if ( metadata_type=="all"){
    slots = c("mean_mat","var_mat","ratio_mat","sscore")
    slots = slots[slots %in% names(get_meta_slot(data_obj,"group_info"))] #only existing slots
  }else{
    slots = metadata_type
  }
  
  group_info = get_meta_slot(data_obj,"group_info")
  for (slot in slots){
    data_mat = group_info[[slot]]
    #filter mapping
    all_mapping = gene_mapping_table %>% 
      dplyr::select(all_of(c(from,to))) %>% 
      dplyr::distinct() %>% 
      dplyr::filter(.[[from]] %in% colnames(data_mat)) %>% #filter by feature names
      tidyr::drop_na(all_of(to)) #drop features without mapping
    #transform to mapping in the tibble
    data_mat = data_mat[,match(all_mapping[[from]],colnames(data_mat))]
    #transform mapping 
    fac_mat = Matrix::fac2sparse(factor(all_mapping[[to]], levels = unique(all_mapping[[to]]))) 
    if(multi_mapping=="mean"){
      fac_mat = sweep_sparse(fac_mat, margin=1, stats = Matrix::rowSums(fac_mat), fun = "/")
    }
    #data_mat = (data_mat %*% Matrix::t(fac_mat)) %>% 
    #  magrittr::set_colnames(unique(all_mapping[[to]]))
    #get_meta_slot(data_obj,"group_info")[[slot]] = data_mat
    group_info[[slot]] = (data_mat %*% Matrix::t(fac_mat)) %>% 
      magrittr::set_colnames(unique(all_mapping[[to]])) 
  }
  data_obj = data_obj %>% set_meta_slot(slot="group_info",group_info)
  return(data_obj)
}
