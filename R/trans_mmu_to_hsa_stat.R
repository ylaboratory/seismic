#' Simply translate the gene names of the  primary information or specificity score calculated via `cal_stat()` (same as the original gene names in the expression matrix) 
#' to human (or other format of) gene id for downstream aim.
#' Mapping translation can be also done between symbol, ensembl id, etc
#' The gene mapping table is also provided as a data frame named `mmu_hsa_mapping`. In this data frame only genes with entrez ID are kept. You can also make your own mapping table.
#' Features with no mapping will be discarded.
#' Because calculation of the variance metric will require extra calculation based on the original data matrix, the translated variance matrix will not be generated.
#' 
#'
#' @param data_obj SingleCellExperiment object
#' @param gene_mapping_table A data frame or tibble to represent gene mapping. 
#' @param from The column (as in gene_mapping_table) to translate the feature names from. 
#' @param to The column (as in gene_mapping_table) to translate the feature names to.
#' @param multi_mapping How to translate metrics when one feature is mapped to several? By "mean" value or "sum" value? (mean for default)
#' @param seismic_data_type The seismic data slots (metadata(data_obj)[["seismic.data"]]) you are translating. By default all seismic analysis data (cell type-level gene mean, ratio of expression and specificity score) will be translate. 
#' You can also specify any one ( or several) if you only would like to translate it (them). The cell-type level gene varaiance 
#' @return A SingleCellExperiment object or a tibble. With seismic analysis data slots are translated into the desired gene names and previous information is kept in the "original_group_info"
#' @export
#' 
trans_mmu_to_hsa_stat = function(data_obj, gene_mapping_table, from, to, multi_mapping = "mean", seismic_data_type = "all") {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  #check column names 
  if (!all(c(from,to) %in% colnames(gene_mapping_table))) {
    stop("The gene_mapping_table is not fine or the columns of 'from' and 'to' do not exist in the table")
  }
  #check which seismic data slot should be used
  #if(!meta_slot_is_null(data_obj, "original_group_info")){
    #group_info = get_meta_slot(data_obj, "original_group_info")
  if(!seismic_slot_is_null(data_obj, "original_group_info")){
    group_info = get_seismic_slot(data_obj, "original_group_info")
  }else{
    #if (meta_slot_is_null(data_obj, "group_info")){
    if (seismic_slot_is_null(data_obj, "group_info")){
      stop("You should run cal_stat() and(or) cal_sscore() first")
    }
    #group_info = get_meta_slot(data_obj, "group_info")
    #data_obj = data_obj %>% set_meta_slot(slot="original_group_info",group_info)  #set as another slot
    group_info = get_seismic_slot(data_obj, "group_info")
    data_obj = data_obj %>% set_seismic_slot(slot="original_group_info",group_info)  #set as another slot
  }
  #check metadata slots 
  if (seismic_data_type!="all" & (!seismic_data_type %in% names(group_info))){
    stop("The seismic data slot does not exist")
  }
  #check feature names
  if ( seismic_data_type =="all"){ 
    #slots = c("mean_mat","ratio_mat","sscore") %>%  .[. %in%  names(get_meta_slot(data_obj,"group_info"))] #only existing slots
    slots = c("mean_mat","ratio_mat","sscore") %>%  .[. %in%  names(get_seismic_slot(data_obj,"group_info"))] #only existing slots
  }else{
    slots = seismic_data_type 
  }
  
  #all_names = slots %>% purrr::map(~colnames(get_meta_slot(data_obj,"group_info")[[.x]])) %>% purrr::reduce(~unique(c(.x,.y))) #all names of genes
  all_names = slots %>% purrr::map(~colnames(get_seismic_slot(data_obj,"group_info")[[.x]])) %>% purrr::reduce(~unique(c(.x,.y))) #all names of genes
  
  if (!any(all_names %in% (gene_mapping_table %>% dplyr::pull(from)))){
    stop("Something's wrong. The feature names don't match the column in the gene_mapping_table")
  }
  
  #checek multi_mapping parameters
  if (!multi_mapping %in% c("mean","sum")) {
    stop("multi_mapping can only be either 'mean' or 'sum'")
  }
  
  #generate a factor matrix for downstream 
  if (any(unlist(map(slots, ~(ncol(group_info[[.x]])!=length(all_names) | any(all_names != colnames(group_info[[.x]]))))))){
    stop("It seems like there's some disconcordance across the symbols of the cell-type level gene information matrix")
  }
  
  #filter mapping
  all_mapping = gene_mapping_table %>% 
    dplyr::select(all_of(c(from,to))) %>% 
    dplyr::distinct() %>% 
    dplyr::filter(.[[from]] %in% all_names) %>% #filter by feature names
    tidyr::drop_na(all_of(to)) #drop features without mapping
  
  #transform mapping 
  gene_fac_mat = Matrix::fac2sparse(factor(all_mapping[[to]], levels = unique(all_mapping[[to]]))) 
  if(multi_mapping=="mean"){
    gene_fac_mat = sweep_sparse(gene_fac_mat, margin=1, stats = Matrix::rowSums(gene_fac_mat), fun = "/")
  }
  
  for (slot in slots){
    data_mat = group_info[[slot]]
    #transform to mapping in the tibble
    data_mat = data_mat[,match(all_mapping[[from]],colnames(data_mat))]
    group_info[[slot]] = (data_mat %*% Matrix::t(gene_fac_mat)) %>% 
      magrittr::set_colnames(unique(all_mapping[[to]])) 
  }
  #obj_log_list = get_meta_slot(data_obj,"obj_log")
  obj_log_list = get_seismic_slot(data_obj,"obj_log")
  obj_log_list[["gene_trans_info"]] = c(from, to)
  group_info = group_info[c("cell_num",slots)]
  data_obj = data_obj %>% 
    #set_meta_slot(slot="group_info",group_info) %>%
    #set_meta_slot(slot="obj_log",obj_log_list)
    set_seismic_slot(slot="group_info",group_info) %>%
    set_seismic_slot(slot="obj_log",obj_log_list)
  
  return(data_obj)
}
