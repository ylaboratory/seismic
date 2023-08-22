#' Plot enrichment results. Return a ggplot object.
#' 
#' @param data_obj SingleCellExperiment object that contains 
#' @param trait Which association result you would like to plot
#' @param show_value Show FDR or Pvalue to plot as y-axis?
#' @param plot_top_option: How to choose cell types to plot? This would be a string parameter passed to dplyr::filter for the data frame containing 4 columns: cell_type, Pvalue, FDR, rank. Also you can filter by cell_type_anno
#' @param group Which annotation column used for grouping?  
#' @param significance To label significance level or not? By default the FDR 0.1 is asterisk, FDR 0.05 is double asterisk, FDR 0.01 is triple asterisk
#' @return A ggplot object
#' @export
#' 
plot_asso = function(data_obj, trait,show_value = "FDR", plot_top_option = NULL ,group = NULL, significance = TRUE  ) {
  if( !show_value %in% c("FDR","Pvalue")){
    stop("Something's wrong with the plot_top_option. Not all columns exist. ")
  }
  if(is.null(S4Vectors::metadata(data_obj)[["association"]][[trait]])){
    stop("The trait-cell type association does not exist")
  }
  if( !is.null(plot_top_option) & !all(all.vars(rlang::parse_expr(plot_top_option)) %in% c("cell_type","Pvalue","FDR","rank",colnames(S4Vectors::metadata(data_obj)[["cell_type_anno"]])))){
    stop("Something's wrong with the plot_top_option. Not all columns exist. ")
  }
  asso_df = S4Vectors::metadata(data_obj)[["association"]][[trait]] %>% 
    dplyr::arrange(Pvalue) %>%
    dplyr::mutate(rank = 1:n())
  #check group
  if (!is.null(group)){
    if(is.null(S4Vectors::metadata(data_obj)[["cell_type_anno"]]) | !group %in%  colnames(S4Vectors::metadata(data_obj)[["cell_type_anno"]])){
      stop("Something's wrong with the group column. To use this you may first do add_ct_anno()")
    }else{
      asso_df = asso_df %>% dplyr::left_join(S4Vectors::metadata(data_obj)[["cell_type_anno"]], by="cell_type")
    }
  }
  
  if(!is.null(plot_top_option)){
    asso_df = asso_df %>% 
      dplyr::filter(!!rlang::parse_expr(plot_top_option))  
  } #filter by option
  
  asso_df = asso_df %>% 
    dplyr::mutate(cell_type = factor(cell_type, levels = rev(cell_type))) %>% #add factor 
    dplyr::mutate(across(all_of(c("Pvalue","FDR")), ~-log10(.)) ) 
    
  #add significance asterisk
  label_less_sig = asso_df %>% dplyr::filter(FDR>= -log10(0.1) & FDR< (-log10(0.05))) %>% dplyr::mutate(across(all_of(c("Pvalue","FDR")), ~.+0.5))
  label_sig = asso_df %>% dplyr::filter(FDR>=-log10(0.05) & FDR< (-log10(0.01)))%>% dplyr::mutate(across(all_of(c("Pvalue","FDR")), ~.+0.5))
  label_most_sig = asso_df %>% dplyr::filter(FDR>=-log10(0.01)) %>% dplyr::mutate(across(all_of(c("Pvalue","FDR")), ~.+0.5))
  
  #plot
  if(is.null(group)){
    plot_obj = ggplot2::ggplot(asso_df, ggplot2::aes(.data[["cell_type"]], .data[[show_value]]))
  }else{
    plot_obj = ggplot2::ggplot(asso_df, ggplot2::aes(.data[["cell_type"]], .data[[show_value]], fill = .data[[group]]))
  }
  plot_obj = plot_obj + 
    ggplot2::geom_bar(stat = "identity") + 
    ggplot2::coord_flip() + 
    ggplot2::geom_text(data = label_less_sig, label = "*") +
    ggplot2::geom_text(data = label_sig, label = "**") +
    ggplot2::geom_text(data = label_most_sig, label = "***") +
    ggplot2::ylab(paste0("-log10(", show_value, ")")) +
    ggplot2::theme_minimal() 
  
  return(plot_obj)
}



