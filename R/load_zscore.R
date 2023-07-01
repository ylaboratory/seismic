#' Load a z-score file from MAGMA output. The output will be saved as a two-column data frame

#' @param file_path Path for the MAGMA .genes.out files
#' @return A data frame with hsa_entrez and z_stat as the two columns (indicating human gene Entrez ID and its respective MAGMA z-score).
#' @import magrittr
#' @export
#' 
load_zscore = function(file_path, name = NULL) {
  if (!file.exists(file_path)){
    stop("The output file does not exist!")
  }
  if (!grepl(pattern = "\\.genes.out",x=file_path)){
    message("The input file doesn't seem to be a gene z-score file from MAGMA, please check it carefully")
  }
  
  if (is.null(name)){
    name = basename(file_path) %>% 
      strsplit(split=".",fixed = T) %>%
      unlist %>% .[1] 
  }
  
  z_score = read.table(file_path, header=T) %>% 
    dplyr::select(GENE, ZSTAT) %>%
    magrittr::set_colnames(c("hsa_entrez",paste0(name,"_zstat")))
  
  return(z_score)
}