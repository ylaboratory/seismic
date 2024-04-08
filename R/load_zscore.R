#' Load z-score file(s) from MAGMA output. 

#' @description Load z-score file(s) from MAGMA output (with .genes.out tails). The file path could be the path of a single file or a directory or even a list of the past file. 
#' For multiple  z-score files, the gene z-score will be combined together in a single data frame.
#' @param file_path Path for the MAGMA .genes.out file (or directory or a list/vector of paths)
#' @return A data frame with hsa_entrez and z_stat as the two columns (indicating human gene Entrez ID and its respective MAGMA z-score).R
#' @export
#' 
load_zscore = function(file_path, name = NULL) {
  if (length(file_path)==1 && file.info(file_path)$isdir){
    file_path = list.files(file_path, full.names=T)
  } #expand a directory
  
  file_path = file_path[file.exists(file_path)] #only keep files that exist
  if (length(file_path)==0){
    stop("The zscore file/directory does not exist!")
  }
  
  file_path = file_path[grepl(pattern = "\\.genes.out",x=file_path)] #only keep files with the right tail 
  if (length(file_path)==0){
    stop("The z-score files are not right specified! They should contain .genes.out tails!")
  }
  
  if (is.null(name)){
    name = basename(file_path) %>% 
      strsplit(split=".",fixed = T) %>%
      purrr::map(~.x[1]) %>%
      unlist 
  }else{
    if (length(name)!=length(file_path)){
      stop("The name argument specified have the different length as the zscore files")
    }
  }
  
  z_score = file_path %>% 
    purrr::map(~read.table(.x, header=T)) %>% 
    purrr::map(~dplyr::select(.x, GENE, ZSTAT)) %>%
    purrr::map2(name, ~magrittr::set_colnames(.x, c("hsa_entrez",paste0(.y,"_zstat")))) %>% 
    purrr::reduce(~full_join(.x,.y,by="hsa_entrez"))
  
  return(z_score)
}


