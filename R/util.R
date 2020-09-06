
#' Calculation of e matrices
#'
#' \code{elk} generates the e matrices.
#'
#' @param n the number of rows
#' @param p the numbwr of columns
#'
#' @return A list of e matrices.
#' \cite{el1} is the matrix of el1 in the paper, which is constructed by column matrix ei. ei is an n-dimensional vector with each component being 0 but its i-th component being 1.
#'
#' \cite{el2} is the matrix of el2 in the paper, which is constructed by column matrix ei. ei is an n-dimensional vector with each component being 0 but its i-th component being 1.
#'
#' \cite{ek1} is the matrix of ek1 in the paper, which is constructed by column matrix ei. ei is an p-dimensional vector with each component being 0 but its i-th component being 1.
#'
#' \cite{ek2} is the matrix of ek2 in the paper, which is constructed by column matrix ei. ei is an p-dimensional vector with each component being 0 but its i-th component being 1.

#' @export
#'
#' @examples
#' elk(4,3)
elk = function(n,p){

  # n,l
  count = 0
  el1 = el2 = matrix(0,n,n*(n-1)/2)

  for(i in (n-1):1){
    temp = matrix(0,n,i)
    temp[n-i,] = 1
    el1[,(count+1):(count+i)] = temp
    el2[(n-i+1):n,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }

  # p,k
  count = 0
  ek1 = ek2 = matrix(0,p,p*(p-1)/2)
  for(i in (p-1):1){
    temp = matrix(0,p,i)
    temp[p-i,] = 1
    ek1[,(count+1):(count+i)] = temp
    ek2[(p-i+1):p,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }
  return(list(el1 = el1, el2 = el2, ek1 = ek1, ek2 = ek2))
}

## Clusterpath preprocessing
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}

vec2tri <- function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}

