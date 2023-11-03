#' Print expression tables for MAGMA-specificity or FUMA model - which can be taken by MAGMA. 
#' 
#' This function will print out expression table for FUMA and MAGMA analysis. Typically the output will including two tables: one main table including the expression gene set (MAGMA) or 
#' expression covriate file (FUMA), which could be taken by MAGMA software for cell type enrichment. The cell type are recoded in to "cluster.1, ..., cluster.n"-like format since MAGMA software
#' may ignore name separator or extreme long string. The auxiliary table contains the mapping from the recoded cell type names to the real cell type names.
#' 
#' @param data_obj SingleCellExperiment object. Note that MAGMA-specificity takes CPM assay and FUMA takes logCPM (or logcounts) as input. 
#' @param group The column name in metadata for grouping
#' @param gene_id_from For the feature names, what kind of gene id is it using? This should be one of: "mmu_symbol", "mmu_ensembl", "mmu_entrez", "hsa_symbol", "hsa_ensembl", "hsa_entrez".
#' @param table_type specify either "MAGMA" or "FUMA" data tables to be printed.
#' @param main_table_path path of the main table.
#' @param aux_table_path Path of the auxiliary table. If this is not specified, not auxiliary tables will be printed out.
#' @param meta_data Metadata dataframe. By default it will be rowData(data_obj)
#' @param assay_name The name of the assay to extract. Recommend CPM for MAGMA-specificity and logCPM for FUMA. If this was not specified, assay_name will be dependent on the table_type parameter.
#' @param filter_rare_ct If cell groups with too few cells should be filtered out. By default this is set to TRUE because these groups may be contaminants or something else. 
#' @param filter_thres If filter_rare_ct is TRUE, what cells we are going to filter out?
#' 
#' @export
#' 
print_exp_tbl = function(data_obj, group, gene_id_from, table_type,main_table_path,  aux_table_path = NULL, meta_data = NULL, assay_name = NULL, filter_rare_ct = TRUE, filter_thres = 20) {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  
  #check table type
  if (!table_type %in% c("MAGMA","FUMA")){
    stop("table_type could be only either 'MAGMA' or 'FUMA'")
  }
  
  #check assay name
  if(table_type=="MAGMA"){
    if (is.null(assay_name)){
      assay_name = "CPM"
    }
    if (assay_name!="CPM"){
      warning(paste0("You are printing out MAGMA gene set tables with ",assay_name," . Better use CPM"))
    }
  }else{
    if (is.null(assay_name)){
      assay_name = "logCPM"
    }
    if (!assay_name %in% c("logCPM","logcounts")){
      warning(paste0("You are printing out FUMA gene set tables with ",assay_name," . Better use logCPM or logcounts"))
    }
  }
  message(paste0("Printing ",table_type, " table with ",assay_name," assay"))
  
  #data("mmu_hsa_mapping")
  #check gene id format
  if(!gene_id_from %in% colnames(mmu_hsa_mapping)){
    stop("The gene_id_from parameter does not match the exisiting mapping table")
  }

  #compute cell type mean
  data_obj = cal_stat(data_obj, meta_data = meta_data, group = group , assay_name = assay_name,mean_only=TRUE)
  if (gene_id_from != "hsa_entrez"){
    data_obj = trans_mmu_to_hsa_stat(data_obj, gene_mapping_table = mmu_hsa_mapping, from = gene_id_from, to = "hsa_entrez")
  }
  mean_mat = metadata(data_obj)[["group_info"]][["mean_mat"]]
  
  #main table directory
  if (!dir.exists(base::dirname(main_table_path))) {
    message(paste0("The directory of the main table path:", base::dirname(main_table_path) ," does not exist... cerated one.!"))
    dir.create(base::dirname(main_table_path), recursive = TRUE)
  }
  
  #if output MAGMA table
  if(table_type=="MAGMA"){
    main_tbl = sweep_sparse(mean_mat, margin=2, stats = colSums(mean_mat),fun="/") %>%
      as.matrix() %>%
      Matrix::t() %>% 
      dplyr::as_tibble(rownames = "hsa_entrez") %>%
      dplyr::pivot_longer(!hsa_entrez, names_to = "cell_type", values_to = "specificity") %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_max(specificity, prop=0.1) %>% 
      dplyr::summarize( genes = paste(hsa_entrez, collapse = " ")) %>%
      dplyr::mutate(cell_type = paste0("cluster.",1:n()))
    write.table(main_tbl,file=main_table_path, col.names = F, row.names = F, sep=" ",quote=F)
  }else{
    #if output FUMA table
    main_tbl = rbind(mean_mat, Matrix::colMeans(mean_mat)) %>%
      magrittr::set_rownames(c(rownames(mean_mat), "Average")) %>%
      as.matrix() %>%
      Matrix::t() %>% 
      dplyr::as_tibble(rownames = "GENE") %>%
      magrittr::set_colnames(c("GENE", paste0("cluster.",1:(ncol(.)-2) ), "Average" ))
    write.table(main_tbl,file=main_table_path, col.names =  T, row.names = F, sep=" ",quote=F)
  }
  message(paste0("The main table has been printed out as ",main_table_path))
  
  
  #if print auxiliary file
  if(!is.null(aux_table_path)){
    aux_tbl = dplyr::tibble(cell_type = rownames(mean_mat)) %>%
      dplyr::mutate(endoced_name = paste0("cluster.",1:n())) 
    #print the table
    if (!dir.exists(base::dirname(aux_table_path))) {
      message(paste0("The directory of the auxiliary table path:", base::dirname(aux_table_path) ," does not exist... cerated one.!"))
              dir.create(base::dirname(aux_table_path), recursive = TRUE)
    }
    write.table(aux_tbl,file=aux_table_path, col.names =  T, row.names = F, sep=" ",quote=F)
    message(paste0("The auxiliary table has been printed out as ",aux_table_path))
  }
}

