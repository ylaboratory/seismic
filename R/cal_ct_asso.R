#' Calculate association between a trait (disease) and a cell type using specificity score after computing specificity score.
#' Genes with too low expression value or only minimum expression in a few cells in the data set should be filter out. 
#' To specify how to filter genes based on these gene annotation, please use a dplyr::filter expression in string.
#' 

#' @param data_obj A SingleCellExperiment object
#' @param gene_zscore_df A data frame containing human entrez id and trait MAGMA z-score. By default they are the first two columns. 
#' @param gene_filter_setting A dplyr::filter expression string for gene filter across the "gene_info" slot in the object meta data. Not setting this indicating that all genes are going to be used.
#' @param asso_model A linear model or spearman correlation ("linear" or "spearman") to test the association between cell type specificity score and MAGMA z-score. 
#' @return A data frame (in metadata(data_obj)[["association"]] slot) containing association p_value, FDR. 
#' @export
#' 
cal_ct_asso = function(data_obj, gene_zscore_df, gene_filter_setting=NULL, asso_model = "linear") {
  #check data
  #if (is.null(S4Vectors::metadata(data_obj)[["group_info"]][["sscore"]]) ){
  if (is.null(get_meta_slot(data_obj,"group_info")[["sscore"]]) ){
    stop("You should first calculate specificity score ")
  }
  #gene filter setting
  if ((!is.null(gene_filter_setting))){
    #if (is.null(S4Vectors::metadata(data_obj)[["gene_info"]])){
    if (meta_slot_is_null(data_obj,"gene_info")){
      stop("Please add gene annotation using 'add_gene_anno()' or 'add_glob_stat()'.")
    #}else if(!all(all.vars(rlang::parse_expr(gene_filter_setting)) %in% colnames(S4Vectors::metadata(data_obj)[["gene_info"]]))){
    }else if(!all(all.vars(rlang::parse_expr(gene_filter_setting)) %in% colnames(get_meta_slot(data_obj,"gene_info")))){
      stop("Something's wrong with the gene_filter_setting. Not all columns exist. ")
    }
  }
  # check model
  if (!asso_model %in% c("linear","spearman")){
    stop("Please indicate the right model.")
  }
  #check if gene_zscore data look good
  #if(length( intersect( colnames(S4Vectors::metadata(data_obj)[["group_info"]][["sscore"]]), gene_zscore_df[[1]]) )==0){
  if(length( intersect( colnames(get_meta_slot(data_obj,"group_info")[["sscore"]]), gene_zscore_df[[1]]) )==0){
    stop("The gene_zscore_df have no genes mapping to the current specificity score gene entry, you may first map genes to the same gene id type.")
  }
  
  #check gene mapping
  if(!any(as.character(gene_zscore_df[[1]]) %in% as.character(colnames(get_meta_slot(data_obj,"group_info")[["sscore"]])))){
    stop("Specificity score do not have the same gene id type as the z_score_df, you may first do mapping between them using trans_mmu_to_hsa_stat()")
  }
  
  #expand to a tibble
  sscore_tb = get_meta_slot(data_obj,"group_info")[["sscore"]] %>%
    as.matrix() %>%
    dplyr::as_tibble(rownames = "cell_type") %>%
    tidyr::pivot_longer(!cell_type,names_to = "gene_name",values_to = "sscore")
  
  #filter genes
  if(!is.null(gene_filter_setting) & !meta_slot_is_null(data_obj,"gene_info")){
    model_gene = get_meta_slot(data_obj,"gene_info") %>%
      dplyr::filter(!!rlang::parse_expr(gene_filter_setting)) %>%
      dplyr::pull(gene_name)
  }else{
    model_gene = colnames(get_meta_slot(data_obj,"group_info")[["sscore"]])
  }
  
  sscore_tb = sscore_tb %>% 
    dplyr::filter(gene_name %in% model_gene) %>%
    dplyr::inner_join(gene_zscore_df %>% dplyr::mutate(hsa_entrez = as.character(hsa_entrez)), by= c("gene_name"="hsa_entrez")) 
  
  #association
  if(asso_model=="linear"){
    asso_res = sscore_tb %>%
      dplyr::group_by(across(all_of("cell_type"))) %>%
      dplyr::summarise( across(colnames(gene_zscore_df)[-1],  ~lm_pvalue(.x, sscore))) 
  }else{
    asso_res = sscore_tb %>%
      dplyr::group_by(across(all_of("cell_type"))) %>%
      dplyr::summarise( across(colnames(gene_zscore_df)[-1],  ~spearman_pvalue(.x, sscore)))
  }
  
  #break to list
  trait_name = gsub(colnames(gene_zscore_df)[-1], pattern="_zstat",replacement = "")
  asso_res = purrr::map(2:ncol(asso_res), ~asso_res[,c(1,.x)])
  asso_res = asso_res %>% purrr::set_names(trait_name) 
  
  #add FDR
  asso_res = asso_res %>%
    purrr::map(~magrittr::set_colnames(.x, c("cell_type","Pvalue"))) %>%
    purrr::map(~dplyr::mutate(.x,FDR = p.adjust(Pvalue,method = "fdr"))) %>%
    purrr::map(~as.data.frame(.x))
  
  #output
  #S4Vectors::metadata(data_obj)[["association"]] = append(S4Vectors::metadata(data_obj)[["association"]],asso_res)
  #only keep non-duplicated traits
  #S4Vectors::metadata(data_obj)[["association"]] = S4Vectors::metadata(data_obj)[["association"]][!duplicated(names( S4Vectors::metadata(data_obj)[["association"]]), fromLast=TRUE)]
  data_obj = add_ct_asso(data_obj, asso_res, names(asso_res), asso_model)
  
  return(data_obj)
}



