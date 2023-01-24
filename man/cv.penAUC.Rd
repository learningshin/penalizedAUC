\name{cv.penAUC}
\alias{cv.penAUC}

\title{ Cross-validation for the penalized AUC-optimizing classification}

\description{
This computes the cross-validated AUC for the lasso- or ridge-penalized AUC-classifier over a grid of values for the tuning parameter, and AUC-optimizing information criterion (AUCIC) for the scad-penalized AUC classifier to choose an optimal tuning parameter.
}
\usage{
cv.penAUC(x, y, nfold = 5, loss = c("exponential", "logistic"), penalty = c("lasso", "ridge", "scad"), criteria = NULL, lambda.set = NULL, learn.rate = 0.01, eps = 1e-04, max.iter = 1000, a = 3.7, seed = 10)
}

\arguments{
  \item{x}{ The design matrix (n x p) with n rows (observations) and p columns (variables). }
  \item{y}{ The \code{{-1, +1}} valued response variable. }
  \item{nfold}{ The number of cross-validation folds. Default: 5. }
  \item{loss}{ Either 'exponential' or 'logistic' loss function. }
  \item{penalty}{ The penalty to be applied to the model. Either 'lasso', 'ridge' or 'scad'. }
  \item{criteria}{ How models will be evaluated. Default is 'CV' for the lasso- or ridge-penalized AUC classifier and 'AUCIC' for the scad-penalized AUC classifier. 'AUCIC' is only applicable to the scad-penalized AUC classifier. }
  \item{lambda.set}{ A vector of lambdas. Default for lambdas is 1.6^(-9:4). }
  \item{learn.rate}{ The value of learning rate, dafault: 1e-2. }
  \item{eps}{ Some small values for the convergence threshold, default: 1e-4. }
  \item{max.iter}{ The maximum value of iterations; default: 1e3.}
  \item{a}{ The tuning parameter for SCAD penalty, default: 3.7. }
  \item{seed}{ The value of seed }
}

\value{
Returns an object of class 'cv.penAUC', a list with the following components:
\itemize{
\item{opt.lambda:} the value of tuning parameter that minimizes the averaged cross-validation risk or AUCIC.
\item{loss:} type of loss function
\item{penalty:} type of penalty function
\item{criteria:} type of criteria for model selection
}
}

\references{
Hyungwoo Kim and Seung Jun Shin, "Variable Selection in AUC-Optimizing Classification"
}

\author{
Hyungwoo Kim and Seung Jun Shin
}

\examples{
library(MASS)

# Data generation
set.seed(3)
n=100; p=10; k=3
class.p  = 0.9
y <- c(rep(1,class.p*n), rep(-1,n-class.p*n))

x <- matrix(NA, n, p)
x[y==1, ] <- mvrnorm(class.p*n, mu=c(rep(0.8,k), rep(0,p-k)), Sigma=diag(1,p))
x[y==-1, ] <- mvrnorm(n-class.p*n, mu=c(rep(-0.8,k), rep(0,p-k)), Sigma=diag(1,p))


############## Find optimal regularization parameter
eps=1e-4; max.iter=1e3; learn.rate=0.01; a=3.7

# Lasso (criteria = 'CV')
lambda.set <- c(1.6^(-8:3)) # grid lambda
cv.lasso.obj <- cv.penAUC(x, y, nfold=5, lambda.set=lambda.set, loss="exponential", penalty="lasso",
                          criteria="CV", learn.rate=learn.rate)
cv.lasso.obj$opt.lambda

# SCAD (criteria = 'AUCIC')
cv.scad.obj <- cv.penAUC(x, y, loss="exponential", penalty="scad", criteria="AUCIC", learn.rate=learn.rate)
cv.scad.obj$opt.lambda
}


