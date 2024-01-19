#' Detect the influential genes in the a trait-cell type association pair using DFBETAS
#' 
#' Can only be used if the specified model is "linear model" and the trait-cell type association has been computed.
#' 

#' @param data_obj A SingleCellExperiment object
#' @param gene_zscore_df A gene zscore data frame. The first column is hsa_entrez and the other columns contain the MAGMA gene Z-score of traits.
#' @param trait_name The name of the trait.
#' @param cell_type The cell type.
#' @param pos_genes If only genes with positive DFBETAS value will be pick as the top influential ones?
#' @return A data frame, containing gene DFBETAS value, cell ttpe gene specificity score, and trait gene MAGMA Z-score
#' @export
#' 

gene_inf_measure = function(data_obj, gene_zscore_df, trait_name, cell_type,  pos_genes = T){
  if(!inherits(data_obj,"SingleCellExperiment")){
    stop("object class fault")
  }
  if(!"linear" %in% names(get_meta_slot(data_obj,"association"))){
    stop("Influential gene detection can only used in a linear setting. This should be done first")
  }
  if(!paste0(trait_name, "_zstat") %in% colnames(gene_zscore_df)){
    stop("The specified traits do not exist in the gene_zscore_df")
  }
  asso_df = get_ct_asso(data_obj, trait_name = trait_name, asso_model = "linear") #also do trait name check
  if(!cell_type %in% asso_df[["cell_type"]]){
    stop("The cell type specified does not exist in the association table.")
  }
  tmp_cell_type = cell_type #avoid misuses 
  fdr_value  = asso_df %>% 
    dplyr::filter(cell_type == tmp_cell_type) %>% 
    pull(FDR)
  if(fdr_value > 0.05){
    warning("The association pair you choose does not seem significant enough (FDR 0.05 threshold here), the results may be unrelieable")
  }
  model_gene = get_meta_slot(data_obj,"association")[["model_gene"]][["linear"]]
  sscore_df = get_meta_slot(data_obj,"group_info")[["sscore"]] %>%
    .[which(rownames(.)==cell_type),] %>%
    as.vector() %>%
    dplyr::as_tibble(rownames = "gene_name") %>% 
    dplyr::mutate(gene_name = as.character(gene_name)) %>% 
    dplyr::filter(gene_name %in% as.character(model_gene)) %>% 
    dplyr::left_join(gene_zscore_df %>% dplyr::select(dplyr::all_of(c("hsa_entrez",paste0(trait_name, "_zstat")))) %>% dplyr::mutate(hsa_entrez = as.character(hsa_entrez)), by=c("gene_name"="hsa_entrez")) %>% 
    magrittr::set_colnames(c("hsa_entrez","specificity_score","trait_z_stat")) %>%
    drop_na()
  #lm
  new_lm = lm(trait_z_stat~specificity_score, data = sscore_df)
  #dfbetas
  dfbeta_df = sscore_df %>% 
    cbind(dfbetas(new_lm)[,2]) %>% 
    magrittr::set_colnames(c("hsa_entrez","specificity_score","trait_z_stat","dfbetas")) 
  if(pos_genes){
    dfbeta_df = dfbeta_df %>% dplyr::mutate(influential = ifelse(dfbetas>2/sqrt(nrow(dfbeta_df)), T, F))
  }else{
    dfbeta_df = dfbeta_df %>% dplyr::mutate(influential = ifelse(abs(dfbetas)>2/sqrt(nrow(dfbeta_df)), T, F))
  }
  return(dfbeta_df)
}