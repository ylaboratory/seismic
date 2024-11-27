#' Similar to the base sweep function, but can be used for sparse matrices as well.
#'
#' Enables fast computations across columns or rows (margin=1 and 2 respectively)
#' Only matrix, dgeMatrix, dgTMatrix, dgCMatrix, or dgRMatrix types are supported.
#' adapted from https://stackoverflow.com/questions/55407656/r-sweep-on-a-sparse-matrix
#'
#' @param x A matrix (matrix, dgeMatrix, dgTMatrix, dgCMatrix, or dgRMatrix).
#' @param margin Across columns or rows (margin=1 and 2 respectively)
#' @param stats The column name in metadata for grouping
#' @param fun function to be applied for each item in stats to x (along margin)
#' @return A matrix

sweep_sparse <- function(x, margin, stats, fun = "*") {
  if ((length(stats) != 1) & (length(stats) != dim(x)[margin])) {
    stop("stats does not match the dimension of the selected margin")
  }
  
  if (class(stats)[1] == "matrix"){
    stats = c(stats)
  }
  
  # to deal with a non-sparse matrix
  if (class(x)[1] %in% c("dgeMatrix", "matrix")) {
    res <- sweep(as.matrix(x), margin, stats, fun)
    if (class(x)[1] == "dgeMatrix"){
      res <- as(res, "dgeMatrix")
    }
    return(res)
  }

  #deal with stats
  if (!class(x)[1] %in% c("dgTMatrix", "dgCMatrix", "dgRMatrix")) {
    stop("Only these sparse matrix types are supported: dgTMatrix, dgCMatrix, dgRMatrix")
  }

  if (class(x)[1] == "dgRMatrix") {
    i <- rep(1:x@Dim[1], diff(x@p)) - 1
  } else {
    i <- x@i
  }

  if (class(x)[1] == "dgCMatrix") {
    j <- rep(1:x@Dim[2], diff(x@p)) - 1
  } else {
    j <- x@j
  }

  if (margin == 1) {
    idx <- i + 1
  } else {
    idx <- j + 1
  }

  if (length(stats) == 1) {
    stats <- rep(stats, x@Dim[margin])
  }

  f <- match.fun(fun)
  x@x <- f(x@x, stats[idx])
  return(x)
}
