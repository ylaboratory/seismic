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
ct_asso = function(data_obj, gene_zscore_df, gene_filter_setting=NULL, asso_model = "linear") {
  #check data
  if (is.null(metadata(data_obj)[["group_info"]][["sscore"]]) ){
    stop("You should first calculate specificity score ")
  }
  #gene filter setting
  if ((!is.null(gene_filter_setting))){
    if (is.null(metadata(data_obj)[["gene_info"]])){
      stop("Please add gene annotation using 'add_gene_anno()' or 'add_glob_stat()'.")
    }else if(!all(all.vars(rlang::parse_expr(gene_filter_setting)) %in% colnames(metadata(data_obj)[["gene_info"]]))){
      stop("Something's wrong with the gene_filter_setting. Not all columns exist. ")
    }
  }
  # check model
  if (!asso_model %in% c("linear","spearman")){
    stop("Please indicate the right model.")
  }
  #check if gene_zscore data look good
  if(length( intersect( colnames(metadata(data_obj)[["group_info"]][["sscore"]]), gene_zscore_df[[1]]) )==0){
    stop("The gene_zscore_df have no genes mapping to the current specificity score gene entry, you may first map genes to the same gene id type.")
  }
  
  #check gene mapping
  if(!any(as.character(gene_zscore_df[[1]]) %in% as.character(colnames(metadata(data_obj)[["group_info"]][["sscore"]])))){
    stop("Specificity score do not have the same gene id type as the z_score_df, you may first do mapping between them using trans_mmu_to_hsa_stat()")
  }
  
  #expand to a tibble
  sscore_tb = metadata(data_obj)[["group_info"]][["sscore"]] %>%
    as.matrix() %>%
    dplyr::as_tibble(rownames = "cell_type") %>%
    tidyr::pivot_longer(!cell_type,names_to = "gene_name",values_to = "sscore")
  
  #filter genes 
  model_gene = metadata(data_obj)[["gene_info"]] %>%
    dplyr::filter(!!rlang::parse_expr(gene_filter_setting)) %>%
    dplyr::pull(gene_name)
  sscore_tb = sscore_tb %>% 
    dplyr::filter(gene_name %in% model_gene) %>%
    dplyr::inner_join(gene_zscore_df %>% dplyr::mutate(hsa_entrez = as.character(hsa_entrez)), by= c("gene_name"="hsa_entrez")) 
  
  #association
  if(asso_model=="linear"){
    asso_res = sscore_tb %>%
      dplyr::group_by(across(all_of("cell_type"))) %>%
      dplyr::summarise( across(colnames(gene_zscore_df)[2],  ~lm_pvalue(.x, sscore))) 
  }else{
    asso_res = sscore_tb %>%
      dplyr::group_by(across(all_of("cell_type"))) %>%
      dplyr::summarise( across(colnames(gene_zscore_df)[2],  ~spearman_pvalue(.x, sscore))) 
  }
  
  asso_res = asso_res %>%
    magrittr::set_colnames(c("cell_type","Pvalue"))  %>%
    mutate(FDR = p.adjust(Pvalue,method = "fdr"))
  
  #output
  name = gsub(colnames(gene_zscore_df)[2],pattern="_zstat",replacement = "")
  metadata(data_obj)[["association"]][[name]] = asso_res
  return(data_obj)
}




