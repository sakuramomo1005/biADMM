#' biC-ADMM: Biclustering Algorithm for Model with Compositional Constraints
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
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#' @param niters Iteraion times
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
#' res3 = biC.ADMM(X, nu1, nu2, nu3, gamma_1, gamma_2,
#'  m, phi, niters = 10, tol = 0.0001, weight.scale = 1, output = 0)
#' dim(res3$A)
biC.ADMM = function(X, nu1, nu2, nu3, gamma_1, gamma_2,
                    m = 5, phi = 0.5,
                    prox = 'l2',
                    niters = 1000, tol = 1e-5,
                    weight.scale = 1, output = 1){

  require(reticulate)
  require(cvxbiclustr)
  require(cvxclustr)
  require(Matrix)
  require(MASS)

  n <- dim(X)[1]; p <- dim(X)[2]

  n2 <- n*(n-1)/2
  p2 <- p*(p-1)/2

  elks <- elk(n,p)
  el1 <- elks$el1
  el2 <- elks$el2
  ek1 <- elks$ek1
  ek2 <- elks$ek2

  k_row <- m
  k_col <- m

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

  w_l <- w_row; u_k <- w_col

  A <- matrix(0,n,p)
  v <- matrix(0,p,n2)
  z <- matrix(0,n,p2)
  lambda_1 <- matrix(0,p,n2)
  lambda_2 <- matrix(0,n,p2)
  lambda_3 <- matrix(0,n,1)

  for(iter in 1:niters){

    A_old <- A; v_old <- v; z_old <- z;
    lambda_1_old <- lambda_1;
    lambda_2_old <- lambda_2
    lambda_3_old <- lambda_3

    # update A

    En <- diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
    Ep <- diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)

    M <- diag(1,n,n) + nu1 * En

    N <- nu2 * Ep + nu3 * matrix(1,p,1) %*% t(matrix(1,p,1))

    s <- matrix(1,n,1) + lambda_3 / nu3

    lv <- lambda_1 + nu1 * v
    lz <- lambda_2 + nu2 * z

    C2 <- (el1-el2) %*% t(lv)
    C3 <- lz %*% t(ek1-ek2)
    C4 <- matrix(rep(s,p),n,p) * nu3

    C <- X +  C2 + C3 + C4

    A <- sylvester(M,t(N),C)

    al1 <- t(A) %*% el1
    al2 <- t(A) %*% el2
    ak1 <- A %*% ek1
    ak2 <- A %*% ek2

    # update vz

    if(prox == 'l1'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      temp1 <- 1 - sigma_1/apply(abs(vtemp),2,sum)
      temp1 <- ifelse(temp1 < 0, 0, temp1)
      temp2 <- matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp
      v <- temp2

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      temp3 <- 1 - sigma_2/apply(abs(ztemp),2,sum)
      temp3 <- ifelse(temp3 < 0, 0, temp3)
      temp4 <- matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp
      z <- temp4

    }else if(prox == 'l2'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      temp1 <- 1 - sigma_1/sqrt(apply(vtemp^2,2,sum))
      temp1 <- ifelse(temp1 < 0, 0, temp1)
      temp2 <- matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp
      v <- temp2

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      temp3 <- 1 - sigma_2/sqrt(apply(ztemp^2,2,sum))
      temp3 <- ifelse(temp3 < 0, 0, temp3)
      temp4 <- matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp
      z <- temp4

    }else if(prox == 'l-inf'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      temp1 <- 1 - sigma_1/apply(abs(vtemp),2,sum)
      temp1 <- ifelse(temp1 < 0, 0, temp1)
      temp2 <- matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp
      v <- vtemp - temp2

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      temp3 <- 1 - sigma_2/apply(abs(ztemp),2,sum)
      temp3 <- ifelse(temp3 < 0, 0, temp3)
      temp4 <- matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp
      z <- ztemp - temp4

    }else{
      print('Error: please specify the norms of the proximal mapping')
      break
    }

    # update lambda
    lambda_1 <- lambda_1 + nu1 * (v - al1 + al2)

    # update lambda 2
    lambda_2 <- lambda_2 + nu2 * (z - ak1 + ak2)

    # update lambda 3
    lambda_3 <- lambda_3 + nu3 * (matrix(1,n,1) - A %*% matrix(1,p,1) )

    if(output == 1){
      print('iter')
      print(iter)
      print(paste('A',mean(abs(A - A_old))))
      print(paste('v',mean(abs(v - v_old))))
      print(paste('z',mean(abs(z -z_old))))
      print(paste('1',mean(abs(lambda_1 - lambda_1_old))))
      print(paste('2',mean(abs(lambda_2 - lambda_2_old))))
      print(paste('3',mean(abs(lambda_3 - lambda_3_old))))
    }

    # whether coverage
    if(mean(abs(A - A_old)) < tol &
       mean(abs(v - v_old)) < tol &
       mean(abs(z - z_old)) < tol &
       mean(abs(lambda_1 - lambda_1_old)) < tol &
       mean(abs(lambda_2 - lambda_2_old)) < tol &
       mean(abs(lambda_3 - lambda_3_old)) < tol){
      return(list(A = A,
                  v = v,
                  z = z,
                  lambda_1 = lambda_1,
                  lambda_2 = lambda_2,
                  lambda_3 = lambda_3,
                  niters = iter))
      break
    }
  }

  if(iter == niters){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A,
                v = v,
                z = z,
                lambda_1 = lambda_1,
                lambda_2 = lambda_2,
                lambda_3 = lambda_3,
                niters = iter))
  }

}


