#' Print expression tables for MAGMA-specificity or FUMA model - which can be taken by MAGMA. 
#' 
#' This function will print out expression table for FUMA and MAGMA analysis. Typically the output will including two tables: one main table including the expression gene set (MAGMA) or 
#' expression covriate file (FUMA), which could be taken by MAGMA software for cell type enrichment. The cell type are recoded in to "cluster.1, ..., cluster.n"-like format since MAGMA software
#' may ignore name separator or extreme long string (the order is the same to the cell type presents in the "cell_num" slot of the "group_info" metadata). The auxiliary table contains the mapping from the recoded cell type names to the real cell type names.
#' 
#' @param data_obj SingleCellExperiment object. Note that MAGMA-specificity takes CPM assay and FUMA takes logCPM (or logcounts) as input. 
#' @param table_type specify either "MAGMA" or "FUMA" data tables to be printed.
#' @param main_table_path path of the main table.
#' @param aux_table_path Path of the auxiliary table. If this is not specified, not auxiliary tables will be printed out.
#' 
#' @export
#' 
print_exp_tbl = function(data_obj, table_type,main_table_path,  aux_table_path = NULL) {
  #check input type
  if (!inherits(data_obj, "SingleCellExperiment")) {
    stop("Only SingleCellExperiment class input is accepted.")
  }
  
  #check table type
  if (!table_type %in% c("MAGMA","FUMA")){
    stop("table_type could be only either 'MAGMA' or 'FUMA'")
  }
  
  if (is.null(S4Vectors::metadata(data_obj)[["group_info"]]) | is.null(S4Vectors::metadata(data_obj)[["group_info"]][["mean_mat"]])){
    stop("The 'gene_info' slot of the metadata is empty. Run cal_stat() first.")
  }
  
  mean_mat = metadata(data_obj)[["group_info"]][["mean_mat"]]
  
  #main table directory
  if (!dir.exists(base::dirname(main_table_path))) {
    message(paste0("The directory of the main table path:", base::dirname(main_table_path) ," does not exist... cerated one.!"))
    dir.create(base::dirname(main_table_path), recursive = TRUE)
  }
  
  #if output MAGMA table
  if(table_type=="MAGMA"){
    main_tbl = sweep_sparse(mean_mat, margin=2, stats = Matrix::colSums(mean_mat),fun="/") %>% 
      as.matrix() %>%
      Matrix::t() %>% 
      dplyr::as_tibble(rownames = "hsa_entrez") %>%
      tidyr::pivot_longer(!hsa_entrez, names_to = "cell_type", values_to = "specificity") %>%
      dplyr::group_by(cell_type) %>%
      dplyr::slice_max(specificity, prop=0.1) %>% 
      dplyr::summarize( genes = paste(hsa_entrez, collapse = " ")) %>%
      dplyr::arrange(match(cell_type, rownames(mean_mat))) %>%
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
    write.table(aux_tbl,file=aux_table_path, col.names =  T, row.names = F, sep="\t",quote=F)
    message(paste0("The auxiliary table has been printed out as ",aux_table_path))
  }
}

