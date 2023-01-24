
library(MASS)

set.seed(4)
n=100; p=10; k=3
class.p  = 0.9
y <- c(rep(1,class.p*n), rep(-1,n-class.p*n))

x <- matrix(NA, n, p)
x[y==1, ] <- mvrnorm(class.p*n, mu=c(rep(0.8,k), rep(0,p-k)), Sigma=diag(1,p))
x[y==-1, ] <- mvrnorm(n-class.p*n, mu=c(rep(-0.8,k), rep(0,p-k)), Sigma=diag(1,p))


eps=1e-3; max.iter=1e3; learn.rate=0.01; a=3.7

############## Find optimal regularization parameter
# Ridge (criteria = 'CV')
lambda.set <- c(1.6^(-8:3)) # grid lambda
cv.ridge.obj <- cv.penAUC(x, y, nfold=5, lambda.set=lambda.set, loss="exponential",
                          penalty="ridge", criteria="CV",eps=eps)
cv.ridge.opt.lamb <- cv.ridge.obj$opt.lambda

# Lasso (criteria = 'CV')
cv.lasso.obj <- cv.penAUC(x, y, nfold=5, lambda.set=lambda.set, loss="exponential",
                          penalty="lasso", criteria="CV", learn.rate=learn.rate, eps=eps)
cv.lasso.opt.lamb <- cv.lasso.obj$opt.lambda

# SCAD (criteria = "AUCIC")
cv.scad.obj <- cv.penAUC(x, y, loss="exponential", penalty="scad", criteria="AUCIC")
cv.scad.opt.lamb <- cv.scad.obj$opt.lambda



################ Fit penAUC with optimal lambda
# Ridge
fit.ridge.obj <- penAUC(x, y, loss="exponential", penalty="ridge", lambda=cv.ridge.opt.lamb)

# Lasso
fit.lasso.obj <- penAUC(x, y, loss="exponential", penalty="lasso", lambda=cv.lasso.opt.lamb)

# SCAD
fit.scad.obj <- penAUC(x, y, loss="exponential", penalty="scad", lambda=cv.scad.opt.lamb)

# Comparison
comp.result = data.frame(fit.ridge = fit.ridge.obj$beta, fit.lasso = fit.lasso.obj$beta,
                         fit.scad =fit.scad.obj$beta)
comp.result
