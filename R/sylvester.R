#' Solve Sylvester equation AX + XB = C.
#'
#'
#' \code{sylvester} finds the matrix X that obey this equation, given matrices A, B, and C.
#'
#' @param A Input matrix
#' @param B Input matrix
#' @param C Input matrix
#' @param tol The tolerance of the algorithm
#'
#' @return X Output matrix
#' @export
#'
#' @examples
#' set.seed(123)
#' A <- matrix(rnorm(9),3,3)
#' B <- matrix(rnorm(9),3,3)
#' Xtrue <- matrix(rnorm(9),3,3)
#' C <- A %*% Xtrue + Xtrue %*% B
#' X <- sylvester(A,B,C)
#' X - Xtrue

sylvester <- function(A, B, C, tol = 0.0001){

  require(MASS)
  require(Matrix)

  A1 <- Schur(A)
  Q1 <- A1$Q; R1 = A1$T

  A2 <- Schur(B)
  Q2 <- A2$Q; R2 = A2$T
  C <- t(Q1) %*% C %*% Q2

  Rsq <- R1 * R1
  I <- diag(dim(A)[1])

  k <- 1
  n <- dim(R2)[1]
  p <- dim(R1)[2]

  X <- matrix(0,p,n)
  while(k < n + 1){
    if(k < n){
      if(abs(R2[k+1, k]) < tol){
        left <- R1 + R2[k,k] * I
        right <- C[,k]
        temp <- matrix(0, dim(X)[1],1)
        if(k == 1){
          temp <- temp
        }else{
          temp <- (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        temp <- matrix(temp, dim(C)[1],1)
        X[,k] <- ginv(left) %*% (right - temp)
        k <- k+1
      }else{
        r11 <- R2[k,k]
        r12 <- R2[k, k+1]
        r21 <- R2[k+1, k]
        r22 <- R2[k+1, k+1]
        temp2 <- matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 <- temp2
          temp3 <- temp3
        }else{
          temps <- X[,1:(k-1)] %*% matrix(R2[1:(k-1),k:(k+1)],k-1,2)
          temp2 <- temps[,1]
          temp3 <- temps[,2]
        }
        b1 <- C[,k] - temp2
        b2 <- C[,k+1] - temp3
        b1_prime <- R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime <- R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime <- matrix(0, dim(X)[1],2)
        b_prime[,1] <- b1_prime
        b_prime[,2] <- b2_prime
        X[,k:(k+1)] <- ginv(R1 %*% R1 + (r11 + r22) * R1 +
                              (r11*r22 - r12*r21) * I) %*% b_prime
        k <- k+2
      }
    }else{
      if(abs(R2[1, k]) > tol){
        left <- R1 + R2[k,k] * I
        right <- C[,k]
        temp <- matrix(0, dim(X)[1],1)
        if(k == 1){
          temp <- temp
        }else{
          temp <- (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        temp <- matrix(temp, dim(C)[1],1)
        X[,k] <- ginv(left) %*% (right - temp)
        k <- k+1
      }else{
        R22 <- R2
        R22 <- cbind(R2, rep(0,dim(R2)[1]))
        R22 <- rbind(R22,rep(0,dim(R2)[1]+1))
        r11 <- R22[k,k]
        r12 <- R22[k, k+1]
        r21 <- R22[k+1, k]
        r22 <- R22[k+1, k+1]
        temp2 <- matrix(0, dim(X)[1],1)
        temp3 <- matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 <- temp2
          temp3 <- temp3
        }else{
          temps <- X[,1:(k-1)] %*% matrix(R22[1:(k-1),k:(k+1)],k-1,2)
          temp2 <- temps[,1]
          temp3 <- temps[,2]
        }

        b1 <- C[,k] - temp2
        b2 <- - temp3
        b1_prime <- R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime <- R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime <- matrix(0, dim(X)[1],2)
        b_prime[,1] <- b1_prime
        b_prime[,2] <- b2_prime
        GOD <- ginv(R1 %*% R1 + (r11 + r22) * R1 +
                      (r11*r22 - r12*r21) * I) %*% b_prime
        X[,k] <- GOD[,1]
        k <- k+2
      }
    }
  }
  return(Q1 %*% X %*% t(Q2))
}
