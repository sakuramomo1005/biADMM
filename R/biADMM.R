#' Algorithms for Convex Biclustering
#'
#' @param X
#' @param nu1
#' @param nu2
#' @param gamma_1
#' @param gamma_2
#' @param kk
#' @param phi
#' @param niter
#' @param tol
#' @param output
#'
#' @return
#' @export
#'
#' @examples
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80)
#' # set parameters
#' nu1 = nu2 = gamma_1 = gamma_2 = 0.1
#' kk = 5
#' phi = 0.5
#' # biADMM algorithm
#' res1 = biADMM(X, nu1, nu2, gamma_1, gamma_2,
#'  kk, phi, niter = 1000, tol = 0.0001, output = 0)

biADMM = function(X, nu1, nu2,
                  gamma_1, gamma_2,
                  kk=5, phi=0.5,niter = 1000,tol = 0.1,output = 1){

  require(reticulate)
  require(cvxbiclustr)
  require(Matrix)
  require(MASS)

  n=dim(X)[1]; p=dim(X)[2]

  n2 = n*(n-1)/2
  p2 = p*(p-1)/2

  elks = elk(n,p)
  el1 = elks$el1
  el2 = elks$el2
  ek1 = elks$ek1
  ek2 = elks$ek2

  k_row = k_col = kk

  w_row <- kernel_weights(t(X), phi/p)
  w_col <- kernel_weights(X, phi/n)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)

  w_l = w_row; u_k = w_col

  A = matrix(0,n,p)
  v = matrix(0,p,n2)
  z = matrix(0,n,p2)
  lambda_1 = matrix(0,p,n2)
  lambda_2 = matrix(0,n,p2)

  for(iter in 1: niter){

    A_old = A; v_old = v; z_old = z; lambda_1_old = lambda_1; lambda_2_old = lambda_2

    # update A

    En = diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
    Ep = diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)

    M = diag(1,n,n) + nu1 * En

    N = nu2 * Ep

    lv = lambda_1+ nu1 * v
    lz = lambda_2 + nu2 * z

    C2 = (el1-el2) %*% t(lv)
    C3 = lz %*% t(ek1-ek2)
    C = X +  C2 + C3

    A = sylvester(M,t(N),C)

    al1 = t(A) %*% el1; al2 = t(A) %*% el2
    ak1 = A %*% ek1; ak2 = A %*% ek2


    # update vz
    sigma_1 = gamma_1 * w_l/nu1
    vtemp = al1 - al2 - 1/nu1 * lambda_1

    temp1 = ifelse((1 - sigma_1/sqrt(apply(vtemp^2,2,sum))) < 0, 0,1 - sigma_1/sqrt(apply(vtemp^2,2,sum)))
    temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp
    v = temp2

    ztemp = ak1 - ak2 - 1/nu2 * lambda_2
    sigma_2 = gamma_2 * u_k/nu2

    temp3 = ifelse((1 - sigma_2/sqrt(apply(ztemp^2,2,sum))) < 0, 0,1 - sigma_2/sqrt(apply(ztemp^2,2,sum)))
    temp4 = matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp

    z = temp4

    # update lambda
    lambda_1 = lambda_1 + nu1 * (v - al1 + al2)

    # update lambda 2
    lambda_2 = lambda_2 + nu2 * (z - ak1 + ak2)
    if(output == 1){
      print('iter')
      print(iter)

      print(paste('A',mean(abs(A - A_old))))
      print(paste('v',mean(abs(v - v_old))))
      print(paste('z',mean(abs(z -z_old))))
      print(paste('l',mean(abs(lambda_1 - lambda_1_old))))
      print(paste('2',mean(abs(lambda_2 - lambda_2_old))))
    }


    # whether coverage
    if(mean(abs(A - A_old)) < tol &
       mean(abs(v - v_old)) < tol&
       mean(abs(z - z_old)) < tol &
       mean(abs(lambda_1 - lambda_1_old)) < tol &
       mean(abs(lambda_2 - lambda_2_old)) <tol){
      return(list(A = A, v = v, z = z,
                  lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
      break
    }
  }

  if(iter == niter){
    #print(paste('not converge within',iter, 'times'))
    return(list(A = A, v = v, z = z,
                lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
  }
}
