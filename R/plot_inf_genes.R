#' Plot influential genes for a given trait and cell type after running find_inf_genes().
#'
#' @param inf_df A data.frame or data.table of influential gene scores output
#' by seismic. Must contain columns that correspond to genes, seismic
#' specificity scores, MAGMA trait z-scores, dfbeta values, and a Boolean
#' column indicating influential genes.
#' @param num_labels Number of top influential genes to label (can also specify "all") with name
#' in provided 'gene_col'. When given specific numbers, the top influential genes (by dfbetas)
#' will be labeled. Only genes that have `is_influential==T` will be labeled.
#' @param consider_neg_zstat A Boolean value indicating whether genes with negative
#' z-scores should be considered (aka for coloring and gene labels). Defaults to FALSE.
#' @param consider_neg_dfbetas A Boolean value indicating whether genes with negative
#' dfbetas should be considered (aka for coloring and gene labels). Defaults to FALSE.
#' @param repel A Boolean value indicating whether or not to use the ggrepel
#' package for clearer labeling. Defaults to TRUE.
#' @param gene_col A character string containing the name of the gene identifier
#' column. Defaults to 'gene' the column output by find_inf_genes().
#' @param spec_col A character string containing the name of the specificity
#' column. Defaults to 'specificity' the column output by find_inf_genes().
#' @param trait_col A character string containing the name of the MAGMA z-score
#' column. Defaults to 'zstat' the column output by find_inf_genes().
#' @param df_betas_col A character string containing the name of the dfbetas
#' column. Defaults to 'dfbetas' the column output by find_inf_genes().
#' @param indicator_col A character string containing the name of the influential
#' gene indicator column.
#' Defaults to 'is_influential' the column output by find_inf_genes().
#' @return A ggplot object of geom_point.
#' @export
#'
plot_inf_genes <- function(inf_df, num_labels = 10,
                           consider_neg_zstat = F, consider_neg_dfbetas = F, repel = T,
                           gene_col = "gene", spec_col = "specificity",
                           trait_col = "zstat", df_betas_col = "dfbetas",
                           indicator_col = "is_influential") {
  lab = zstat = is_influential = dfbetas = plot_gene = rownum = NULL # due to non-standard evaluation notes in R CMD check

  # check that the influential gene matrix is valid
  for (col in c(gene_col, spec_col, trait_col, df_betas_col, indicator_col))
  {
    if (!(col %in% names(inf_df))) {
      stop("Missing ", col, " from the influential gene input. Please correct
           the influential data frame format and try again.")
    }
  }

  inf_df <- as.data.table(inf_df)
  setnames(
    inf_df,
    c(gene_col, spec_col, trait_col, df_betas_col, indicator_col),
    c("plot_gene", "specificity", "zstat", "dfbetas", "is_influential")
  )
  inf_df[, lab := ""]

  if (!(consider_neg_zstat)) {
    inf_df[zstat < 0, is_influential := F]
  }

  if (!(consider_neg_dfbetas)) {
    inf_df[dfbetas < 0, is_influential := F]
    # do not want to label genes with negative dfbetas
    inf_df[dfbetas < 0, dfbetas := 0]
  }

  if (num_labels == "all") {
    inf_df[, lab := ifelse(inf_df$is_influential == T, plot_gene, lab)]
  } else if (is.numeric(num_labels)) {
    inf_df <- inf_df[order(-abs(dfbetas)), rownum := 1:.N]
    inf_df[is_influential == F, rownum := NA]
    inf_df[, lab := ifelse(!is.na(rownum) & rownum <= num_labels, plot_gene, lab)]
  } else {
    stop("Only 'all' or numeric values possible for num_labels, not ", num_labels)
  }

  gg <- ggplot2::ggplot(inf_df, ggplot2::aes_string(
    x = "specificity", y = "zstat",
    color = "is_influential", alpha = 0.8,
    label = "lab"
  )) +
    ggplot2::geom_point() +
    ggplot2::ylab("gene risk (z-score)") +
    ggplot2::xlab("specificity score") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("grey", "red"))

  if (repel) {
    gg <- gg + ggrepel::geom_text_repel(
      force = 10,
      alpha = 1,
      segment.square = F,
      segment.inflect = T,
      size = 3,
      max.overlaps = Inf,
      segment.size = 0.1,
      segment.color = "grey10",
      colour = "black"
    )
  } else {
    gg <- gg + ggplot2::geom_text(nudge_x = 0, nudge_y = 0.5, color = "black")
  }

  return(gg)
}
