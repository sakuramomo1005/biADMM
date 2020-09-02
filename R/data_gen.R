#' Simulated data set generation
#'
#' \code{data_gen} generates datasets in the simulation study.
#'
#' @param seed.cluster seed.cluster is used to control the cluster assignment
#' @param seed.data seed.data is used to control the data generation given clustering structure
#' @param n number of subjects
#' @param p number of features
#' @param theta the variance of the x, which is sampled from a normal distribution
#' @param row_group number of row clusters
#' @param col_group number of column clusters
#'
#' @return Output is the simulated dataset
#' @export
#'
#' @examples
#' set.seed(123)
#' data = data_gen(n = 100, p = 80)

data_gen <- function(seed.cluster = 123, seed.data = 654, n, p, theta = 2, row_group = 4, col_group = 4){

  mu <- seq(-10,10,1)
  x <- matrix(0, n, p)

  set.seed(seed.cluster)
  mu_rc <- sample(mu, size = row_group*col_group, replace = T)
  dim(mu_rc) <- c(row_group,col_group )

  row_assign <- c()
  col_assign <- c()
  for(i in 1:n){
    row_assign <- c(row_assign, sample(1:row_group)[1])
  }
  for(i in 1:p){
    col_assign <- c(col_assign, sample(1:col_group)[1])
  }

  ################
  set.seed(seed.data)
  for(i in 1:n){
    for(j in 1:p){
      r <- row_assign[i]
      c <- col_assign[j]
      mu_temp <- mu_rc[r,c]
      x[i,j] <- rnorm(1,mu_temp,theta)
    }
  }

  rownames(x) <- paste('row', row_assign, sep = '')
  colnames(x) <- paste('col', col_assign, sep = '')

  return(x)
}
