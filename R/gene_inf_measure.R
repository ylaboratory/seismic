#' Detect the influential genes in the a trait-cell type association pair using DFBETAS
#' 
#' Can only be used if the specified model is "linear model" and the trait-cell type association has been computed.
#' 

#' @param data_obj A SingleCellExperiment object
#' @param gene_zscore_df A gene zscore data frame. The first column is hsa_entrez and the other columns contain the MAGMA gene Z-score of traits.
#' @param cell_type The cell type.
#' @param trait_name The name of the trait. By default it will be the first trait of the gene_zsore_df.
#' @param pos_genes If only genes with positive DFBETAS value will be pick as the top influential ones?
#' @return A data frame, containing gene DFBETAS value, cell ttpe gene specificity score, and trait gene MAGMA Z-score
#' @export
#' 

gene_inf_measure = function(data_obj, gene_zscore_df,  cell_type, trait_name=NULL,  pos_genes = T){
  if(!inherits(data_obj,"SingleCellExperiment")){
    stop("object class fault")
  }
  #if(!"linear" %in% names(get_meta_slot(data_obj,"association")) | get_meta_slot(data_obj,"obj_log")[["progress"]]!="cal_ct_asso()"){
  if(!"linear" %in% names(get_seismic_slot(data_obj,"association")) | get_seismic_slot(data_obj,"obj_log")[["progress"]]!="cal_ct_asso()"){
    stop("Influential gene detection can only used in a linear model setting. You should first run the cal_ct_asso() function first.")
  }
  if(is.null(trait_name) ){
    trait_name = gsub(colnames(gene_zscore_df)[2],pattern = "_zstat",replacement = "")
  }else{
    if(!paste0(trait_name, "_zstat") %in% colnames(gene_zscore_df)){
      stop("The specified traits do not exist in the gene_zscore_df")
    }
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
  #model_gene = get_meta_slot(data_obj,"obj_log")[["asso_model"]][["linear"]][["model_genes"]]
  model_gene = get_seismic_slot(data_obj,"obj_log")[["asso_model"]][["linear"]][["model_genes"]]
  #sscore_df = get_meta_slot(data_obj,"group_info")[["sscore"]] %>%
  sscore_df = get_seismic_slot(data_obj,"group_info")[["sscore"]] %>%
    .[which(rownames(.)==cell_type),] %>%
    dplyr::as_tibble(rownames = "gene_name") %>% 
    dplyr::mutate(gene_name = as.character(gene_name)) #
  #warning if the number of intersected genes seemed to be wrong
  if(length(intersect(as.character(sscore_df$gene_name), as.character(model_gene)))<length(model_gene)/2){
    warning("It seemed like the genes in the current specificity score matrix do not match that in the previous study.")
  }
    
  sscore_df = sscore_df %>%
    dplyr::filter(gene_name %in% as.character(model_gene)) %>% 
    dplyr::left_join(gene_zscore_df %>% dplyr::select(dplyr::all_of(c("hsa_entrez",paste0(trait_name, "_zstat")))) %>% dplyr::mutate(hsa_entrez = as.character(hsa_entrez)), by=c("gene_name"="hsa_entrez")) %>% 
    magrittr::set_colnames(c("hsa_entrez","specificity_score","trait_z_stat")) %>%
    drop_na()
  #lm
  new_lm = lm(trait_z_stat~specificity_score, data = sscore_df)
  #dfbetas
  dfbeta_df = sscore_df %>% 
    cbind(dfbetas(new_lm)[,2]) %>% 
    magrittr::set_colnames(c("hsa_entrez","specificity_score",paste0(trait_name,"_z_stat"),"dfbetas")) 
  if(pos_genes){
    dfbeta_df = dfbeta_df %>% dplyr::mutate(influential = ifelse(dfbetas>2/sqrt(nrow(dfbeta_df)), T, F))
  }else{
    dfbeta_df = dfbeta_df %>% dplyr::mutate(influential = ifelse(abs(dfbetas)>2/sqrt(nrow(dfbeta_df)), T, F))
  }
  return(dfbeta_df)
}