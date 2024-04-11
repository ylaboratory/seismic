#'sweep_sparse function do the same as the base sweep function, but is adapted to sparse matrix.
#'Fast computation across columns or rows (margin=1 and 2 respectively) for the matrix.
#'Only dgCMatrix and dgRMatrix are accepted.
#' @param x A sparse matrix.
#' @param margin Across columns or rows (margin=1 and 2 respectively)
#' @param fun The column name in metadata for grouping
#' @return A sparse matrix
#' @keywords internal
#' @noRd
sweep_sparse = function(x, margin, stats, fun = "*") {
  #handle error 
  if((length(stats)!=1) & (length(stats)!=dim(x)[margin])){
    stop("stats doesn't match the dimension of the certain margin")
  }
  #to deal with a non-sparse matrix
  if (class(x)[1] %in% c("dgeMatrix","matrix")){
    return(sweep(x,margin,stats,fun))
  }
  if (!class(x)[1] %in% c("dgTMatrix","dgCMatrix","dgRMatrix")){
    stop("Only several saprse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }
  if (class(x)[1]=="dgRMatrix"){
    i = rep(1:x@Dim[1], diff(x@p)) -1
  }else{
    i = x@i
  }
  if (class(x)[1]=="dgCMatrix"){
    j = rep(1:x@Dim[2], diff(x@p)) -1
  }else{
    j = x@j
  }
  if (margin == 1) {
    idx = i + 1
  } else {
    idx = j + 1
  }
  if(length(stats )==1){
    stats = rep(stats, x@Dim[margin])
  }
  f = match.fun(fun)
  x@x = f(x@x, stats[idx])
  return(x)
}


#'This function transform sparsematrix (using such as log functions)
#'input: sparse matrix, margin, stats and function
#'output:sparse matrix.
#' @param x A sparse matrix
#' @param fun A function that only takes one parameter
#' @keywords internal
#' @noRd
transform_sparse = function(x, fun = "log2") {
  if (!class(x)[1] %in% c("dgTMatrix","dgCMatrix","dgRMatrix")){
    stop("Only several saprse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }
  f = match.fun(fun)
  x@x = f(x@x)
  return(x)
}

#'This function use a linear regression to test if two variable are positively correlated 
#'input: sparse matrix, margin, stats and function
#'output:sparse matrix.
#' @param zscore A zscore vector
#' @param metric Value to test association
#' @keywords internal
#' @noRd
lm_pvalue = function(zscore, metric){
  f1 = which(is.finite(metric))
  asso_lm = speedglm::speedlm(zscore[f1]~metric[f1]) #the function will automatically drop the NA values but cannot deal with metric 
  #asso_lm = lm(zscore[f1]~metric[f1])
  lm_sum = summary(asso_lm)$coefficients
  if(lm_sum[2,1]<0) return(1-lm_sum[2,4]/2)
  return(lm_sum[2,4]/2)
}

#'This function use a linear regression to test if two variable are positively correlated 
#'input: sparse matrix, margin, stats and function
#'output:sparse matrix.
#' @param zscore A zscore vector
#' @param metric Value to test association
#' @keywords internal
#' @noRd
spearman_pvalue = function(zscore, metric){
  cor.model = cor.test(z_score, metric, alternative = "greater", method = "spearman")
  return(cor.model$p.value)
}





