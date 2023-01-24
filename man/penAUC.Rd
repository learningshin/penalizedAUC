\name{penAUC}
\alias{penAUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fit the penalized AUC-optimizing classification}

\description{
This provides the fitted coefficients for the AUC-optimizing classification penalized by 'lasso', 'ridge', and 'scad'.
}

\usage{
penAUC(x, y, loss = c("exponential", "logistic"), penalty = c("lasso", "ridge", "scad"), lambda, learn.rate = 0.01, eps = 1e-04, max.iter = 1000, a = 3.7)
}

\arguments{
  \item{x}{ The design matrix (n x p) with n rows (observations) and p columns (variables).}
  \item{y}{ The \code{{-1, +1}} valued response variable. }
  \item{loss}{ Either 'exponential' or 'logistic' loss function. }
  \item{penalty}{ The penalty to be applied to the model. Either 'lasso', 'ridge', or 'scad'. }
  \item{lambda}{ The value of tuning parameter.  }
  \item{learn.rate}{ The value of learning rate, dafault: 1e-2. }
  \item{eps}{ Some small values for the convergence threshold, default: 1e-4. }
  \item{max.iter}{ The maximum value of iterations; default: 1e3. }
  \item{a}{ The tuning parameter for SCAD penalty, default: 3.7. }
}

\value{
Returns an object of class 'penAUC', a list with the following components:
\itemize{
\item{beta:} fitted coefficients
\item{lambda:} value of tuning parameter
\item{loss:} type of loss function
\item{penalty:} type of penalty function
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

set.seed(4)
n=100; p=10; k=3
class.p  = 0.9
y <- c(rep(1,class.p*n), rep(-1,n-class.p*n))
x <- matrix(NA, n, p)
x[y==1, ] <- mvrnorm(class.p*n, mu=c(rep(0.8,k), rep(0,p-k)), Sigma=diag(1,p))
x[y==-1, ] <- mvrnorm(n-class.p*n, mu=c(rep(-0.8,k), rep(0,p-k)), Sigma=diag(1,p))


################ Fit penAUC
fit.lasso.obj <- penAUC(x, y, loss="exponential", penalty="lasso", lambda=.1)
fit.ridge.obj <- penAUC(x, y, loss="exponential", penalty="ridge", lambda=.1)
fit.scad.obj  <- penAUC(x, y, loss="exponential", penalty="scad" , lambda=.1)

# comparison
data.frame(fit.lasso = fit.lasso.obj$beta, fit.ridge = fit.ridge.obj$beta, fit.scad = fit.scad.obj$beta)
}



