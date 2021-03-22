#' biC-ADMM: Biclustering Algorithm for Model with Compositional Constraints (faster version)
#'
#' Same algorithm as \code{biC.ADMM}. Call python code to speed up the running time.
#'
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param nu1 A regularization parameter for row shrinkage
#' @param nu2 A regularization parameter for column shrinkage
#' @param nu3 A regularization parameter for compositional data constrain
#' @param gamma_1 A regularization parameter for row shrinkage
#' @param gamma_2 A regularization parameter for column shrinkage
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param prox The proximal maps. Could calculate L1 norm, L2 norm, or L-infinity, use "l1", "l2", or "l-inf", respectively.
#' @param niters Iteraion times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#' @param weight.scale If weight.scale = 1, the code will make the input data have compositional structure.
#'
#' @return A list of results, containing matrix of A, v, z, lambda1, lambda2, and lambda3
#' @export
#'
#' @examples
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80)
#' # set parameters
#' nu1 = nu2 = nu3 = gamma_1 = gamma_2 = 0.1
#' m = 5
#' phi = 0.5
#' # biADMM algorithm
#' res4 = biC.ADMM.speed(X, nu1, nu2, nu3, gamma_1, gamma_2,
#'  m, phi, niter = 10, tol = 0.0001, weight.scale = 1, output = 0)
#' dim(res4$A)
biC.ADMM.speed = function(X, nu1, nu2, nu3,
                          gamma_1, gamma_2,
                          m, phi,
                          prox = 'l2',
                          niters = 10, tol = 0.1,
                          weight.scale = 1,
                          output = 1){

  require(reticulate)
  require(cvxbiclustr)
  require(cvxclustr)
  require(Matrix)
  require(MASS)

  path <- paste(system.file(package="biclusterADMM"), "biADMM.python.py", sep="/")
  source_python(path)

  n <- dim(X)[1]; p <- dim(X)[2]

  k_row <- m; k_col <- m

  if(weight.scale==1){
    w_row <- kernel_weights(scale(t(X)), phi/p)
    s.X <- scale(X)
    s.X[is.na(s.X)] = 0
    w_col <- kernel_weights(s.X, phi/n) # scale cols to make them comparible
    w_row <- knn_weights(w_row, k_row, n)
    w_col <- knn_weights(w_col, k_col, p)
    w_row <- w_row/sum(w_row)
    w_col <- w_col/sum(w_col)
    w_row <- w_row/sqrt(p)
    w_col <- w_col/sqrt(n)
  } else {
    w_row <- kernel_weights(t(X), phi/p)
    w_col <- kernel_weights(X, phi/n)
    w_row <- knn_weights(w_row, k_row, n)
    w_col <- knn_weights(w_col, k_col, p)
    w_row <- w_row/sum(w_row)
    w_col <- w_col/sum(w_col)
    w_row <- w_row/sqrt(p)
    w_col <- w_col/sqrt(n)
  }

  w_l <- w_row
  u_k <- w_col
  w_l <- matrix(w_l, length(w_l),1)
  u_k <- matrix(u_k, length(u_k),1)

  res <- biADMM_python_compositional(X, nu1, nu2, nu3, gamma_1, gamma_2,
                                     w_l, u_k,
                                     prox,
                                     niters, tol, output = output)

  result <- list(A = res[[1]], v = res[[2]], z = res[[3]],
                 lambda_1 = res[[4]], lambda_2 = res[[5]], lambda_3 = res[[6]],
                 iters = res[[7]])
  return(result)
}
