#' Algorithms for Convex Biclustering
#'
#' @param X
#' @param nu1
#' @param nu2
#' @param gamma_1
#' @param gamma_2
#' @param kk
#' @param phi
#' @param niters
#' @param tol
#' @param output
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
biADMM.speed = function(X,nu1,nu2, gamma_1, gamma_2, kk, phi, niters = 10, tol = 0.1, output = 1){

  path <- paste(system.file(package="biADMM"), "biADMM.python.py", sep="/")
  source_python(path)

  n = dim(X)[1]; p = dim(X)[2]

  k_row <- kk; k_col <-kk
  w_row <- kernel_weights(t(X), phi/p)
  w_col <- kernel_weights(X, phi/n)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)

  w_l = w_row; u_k = w_col
  w_l = matrix(w_l, length(w_l),1)
  u_k = matrix(u_k, length(u_k),1)

  res = biADMM_python(X, nu1, nu2, gamma_1, gamma_2, w_l, u_k, kk, phi, niters, tol, output = 1)

  result = list(A = res[[1]], v = res[[2]], z = res[[3]], lambda_1 = res[[4]], lambda_2 = res[[5]], iters = res[[6]])
  return(result)
}
