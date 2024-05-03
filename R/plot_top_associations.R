#' Plot the top cell type associations for a given trait.
#'
#' @param trait_res A data.frame or data.table object of cell type associations
#' detected by seismic. It is expected that this data object contains a minimum
#' of three columns: cell_type, p_value, and FDR.
#' @param fdr A Boolean value determining if FDR should be plotted. If false
#' p-values are used for the plot. Defaults to printing FDR.
#' @param limit A positive integer limiting the number of top cell types to
#' be shown. Default is 10 cell types.
#' @return A ggplot bar plot object.
#' @export
#'
plot_top_associations <- function(trait_res, fdr = T, limit = 10) {
  score <- FDR <- pvalue <- cell_type <- NULL # due to non-standard evaluation notes in R CMD check

  trait_res <- as.data.table(trait_res)
  xlabel <- ""

  # check for the minimum necessary columns of the
  # trait association matrix
  if (fdr && !("FDR" %in% names(trait_res))) {
    stop("Invalid trait association matrix. Missing FDR column.")
  }
  if (!fdr && !("pvalue" %in% names(trait_res))) {
    stop("Invalid trait association matrix. Missing pvalue column.")
  }
  if (!("cell_type" %in% names(trait_res))) {
    stop("Invalid trait association matrix. Missing cell_type column.")
  }
  if (limit <= 0) {
    stop("Please enter a value greater than 0 for the cell type limit.")
  }

  if (fdr) {
    trait_res[, score := FDR]
    xlabel <- "FDR"
  } else {
    trait_res[, score := pvalue]
    xlabel <- "p-value"
  }

  trait_res[, score := -log10(score)]
  xlabel <- paste0("-log10(", xlabel, ")")

  # filter to the top k values
  trait_res <- trait_res[order(pvalue)] # for ties in fdr
  trait_res <- trait_res[order(-score)]
  trait_res <- utils::head(trait_res, limit)

  trait_res[, cell_type := factor(cell_type, levels = rev(cell_type))]
  gg <- ggplot2::ggplot(trait_res, ggplot2::aes(x = cell_type, y = score)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_classic() +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab(xlabel)

  return(gg)
}
