#' Detect the influential genes in the a trait-cell type association pair using DFBETAS
#' 
#' Can only be used if the specified model is "linear model" and the trait-cell type association has been computed.
#' 

#' @param data_obj A SingleCellExperiment object
#' @param cell_type The cell type.
#' @param trait_name The name of the trait.
#' @param pos_genes If only genes with positive DFBETAS value will be pick as the top influential ones?
#' @return A data frame, containing gene DFBETAS value, cell ttpe gene specificity score, and trait gene MAGMA Z-score
#' @export
#' 

gene_inf_measure = function(data_obj, cell_type, trait_name, pos_genes = T){
  
}