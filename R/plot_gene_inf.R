#' Scatter plot for influential gene testing 
#' 
#' Can only be used if the specified model is "linear model" and the trait-cell type association has been computed.
#' 
#' @param dbetas_df The data frame from previous gene_inf_measure() function. Extra columns can be added to specify the gene annotation in the plot
#' @param gene_label_column The column name to add 
#' @param num_top_gene_label Number of top influential genes (based on DFBETAS values) to label. If it's 0 then no genes are labelled. 
#' @param label_pos_genes If only genes with positive Z score should be labelled? By default this is TRUE.
#' @param extra_labels Extra labels (besides these indicating using "num_top_gene_label" or "label_pos_genes")
#' @param label_repel If the labels will be expelled using the ggrepel. By default this is FALSE. This is recommended only when you have ggrepel installed.
#' @param ... Additional arguments passed on to `geom_text()` or `geom_text_repel()`. This can include any aesthetic properties or parameters supported by `geom_text()`, such as `size`, `color`, `angle`, and so on.
#' @return A ggplot object
#' @export
#' 
plot_gene_inf = function(dfbetas_df, gene_label_column="hsa_entrez", num_top_gene_label = 0, label_pos_genes=T,extra_labels = NULL, label_repel=F, ... ) {
  #check column names
  if(any(!c("hsa_entrez","specificity_score","dfbetas","influential",gene_label_column) %in% colnames(dfbetas_df)) | !any( grepl("_z_stat",colnames(dfbetas_df)) )){
    stop("The data frame does not seem to come from the previous gene_inf_measure() function or the gene_label_column does not exist.")
  }
  trait_column = colnames(dfbetas_df)[which(grepl("_z_stat",colnames(dfbetas_df)))]
  
  if(label_pos_genes){
    dfbetas_df = dfbetas_df %>% dplyr::mutate(influential = ifelse(.data[[trait_column]]>0, influential,FALSE))
  }
  
  if (num_top_gene_label>0){
    dfbetas_df = dfbetas_df %>% 
      dplyr::arrange(-influential, -dfbetas) %>%
      dplyr::mutate(text_label = ifelse(1:n()<=num_top_gene_label, .data[[gene_label_column]],""))
  }
  #deal with extra labels:
  if (is.null(extra_labels)){
    extra_labels = c()
  }else{
    extra_labels = intersect(extra_labels, dfbetas_df[[gene_label_column]])
  }
  if (length(extra_labels)>0){
    dfbetas_df = dfbetas_df %>% 
      dplyr::mutate(text_label = ifelse(.data[[gene_label_column]] %in% extra_labels, .data[[gene_label_column]], text_label)) 
  }
  
 
  #plot
  if("text_label" %in% colnames(dfbetas_df)){
    plot_obj = ggplot2::ggplot(dfbetas_df, ggplot2::aes(x=.data[["specificity_score"]], y = .data[[trait_column]], color=influential,alpha=influential,label=text_label)) 
  }else{
    plot_obj = ggplot2::ggplot(dfbetas_df, ggplot2::aes(x=.data[["specificity_score"]], y = .data[[trait_column]], color=influential,alpha=influential))
  }
    
  
  plot_obj = plot_obj + 
    ggplot2::geom_point() +
    ggplot2::theme_minimal() + 
    ggplot2::scale_color_manual(values =c("grey","red")) +
    ggplot2::scale_alpha_manual(values =c(0.3,0.8))  
  
  if("text_label" %in% colnames(dfbetas_df)){
    plot_obj = plot_obj + geom_point(data = dfbetas_df %>% dplyr::filter(text_label != ""),
                                     aes(x=.data[["specificity_score"]], y = .data[[trait_column]]),
                                     color="black", shape=21, fill=NA,stroke=0.5)
    if(label_repel){
      plot_obj = plot_obj + ggrepel::geom_text_repel( force  = 5,
                                                      color="black",alpha=1,
                                                      segment.square = F,
                                                      segment.inflect =T,
                                                      size=2,
                                                      max.overlaps = Inf,
                                                      ...)
    }else{
      plot_obj = plot_obj + geom_text(color="black",alpha=1,...)
    }
  }
  
  return(plot_obj)
}

