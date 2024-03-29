% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biADMM.compositional.speed.R
\name{biC.ADMM.speed}
\alias{biC.ADMM.speed}
\title{biC-ADMM: Biclustering Algorithm for Model with Compositional Constraints (faster version)}
\usage{
biC.ADMM.speed(
  X,
  nu1,
  nu2,
  nu3,
  gamma_1,
  gamma_2,
  m,
  phi,
  prox = "l2",
  niters = 10,
  tol = 0.1,
  weight.scale = 1,
  output = 1
)
}
\arguments{
\item{X}{The data matrix to be clustered. The rows are the samples, and the columns are the features.}

\item{nu1}{A regularization parameter for row shrinkage}

\item{nu2}{A regularization parameter for column shrinkage}

\item{nu3}{A regularization parameter for compositional data constrain}

\item{gamma_1}{A regularization parameter for row shrinkage}

\item{gamma_2}{A regularization parameter for column shrinkage}

\item{m}{m-nearest-neighbors in the weight function}

\item{phi}{The parameter phi in the weight function}

\item{prox}{The proximal maps. Could calculate L1 norm, L2 norm, or L-infinity, use "l1", "l2", or "l-inf", respectively.}

\item{niters}{Iteraion times}

\item{tol}{Stopping criterion}

\item{weight.scale}{If weight.scale = 1, the code will make the input data have compositional structure.}

\item{output}{When output = 1, print the results at each iteration. No print when output equals other value.}
}
\value{
A list of results, containing matrix of A, v, z, lambda1, lambda2, and lambda3
}
\description{
Same algorithm as \code{biC.ADMM}. Call python code to speed up the running time.
}
\examples{
# generate dataset
set.seed(123)
X = data_gen(n = 100, p = 80)
# set parameters
nu1 = nu2 = nu3 = gamma_1 = gamma_2 = 0.1
m = 5
phi = 0.5
# biADMM algorithm
res4 = biC.ADMM.speed(X, nu1, nu2, nu3, gamma_1, gamma_2,
 m, phi, niter = 10, tol = 0.0001, weight.scale = 1, output = 0)
dim(res4$A)
}
