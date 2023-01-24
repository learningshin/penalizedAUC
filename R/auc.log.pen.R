
auc.log.pen <- function(x, y, penalty=c("lasso","ridge","scad"), lambda,
                        eps=eps, max.iter=max.iter, learn.rate=learn.rate, a=a){

  # scale
  x <- scale(x, scale=F)
  y <- y
  p <- ncol(x)
  i_1 <- which(y==1)
  i_2 <- which(y==-1)
  nm  <- c(length(i_1)*length(i_2))

  # step size, lamb, max.iter
  max.iter <- max.iter
  beta.old = rep(0,p)

  # Fit the AUCpen with exponential loss
  if(penalty == "lasso"){

    for(iter in 1:max.iter){
      t.s <- 1/(10+iter)
      m.ij <- outer( c((y[i_1]*x[i_1,])%*%beta.old), c((y[i_2]*x[i_2,])%*%beta.old), "+")

      x_p = x[i_1,]
      x_m = x[i_2,]
      derv.m <- -(x_p[rep(1:nrow(x_p), times=rep(length(i_2),nrow(x_p))),] - x_m[rep(1:length(i_2), length(i_1)),])
      grad.v <- apply( (exp(-c(t(m.ij)))/(1+exp(-c(t(m.ij))))) * derv.m, 2, mean)

      el.old <- beta.old - t.s* grad.v
      beta.new <- S.f(el.old,lambda,t.s)

      if ( max(beta.new - beta.old) < eps ) break
      beta.old <- beta.new
    }

  } else if(penalty == "ridge"){

    for(iter in 1:max.iter){
      t.s <- 1/(10+iter)
      m.ij <- outer( c((y[i_1]*x[i_1,])%*%beta.old), c((y[i_2]*x[i_2,])%*%beta.old), "+")

      x_p = x[i_1,]
      x_m = x[i_2,]
      derv.m <- -(x_p[rep(1:nrow(x_p), times=rep(length(i_2),nrow(x_p))),] - x_m[rep(1:length(i_2), length(i_1)),])
      grad.v <- apply((exp(-c(t(m.ij)))/(1+exp(-c(t(m.ij))))) * derv.m, 2, mean) + lambda*beta.old
      beta.new <- beta.old - t.s*grad.v

      if ( max(abs(beta.new - beta.old)) < eps ) break
      beta.old <- beta.new
    }

  } else{

    # step size
    t.s <- learn.rate
    beta.c = rep(0,p)
    beta.init <- penAUC(x, y, loss="logistic", penalty="ridge", lambda, learn.rate, eps, max.iter)

    for(i in 1:max.iter){
      beta.old <- beta.init$beta

      for(iter in 1:max.iter){
        t.s <- 1/(iter+10)
        m.ij <- outer( c((y[i_1]*x[i_1,])%*%beta.old), c((y[i_2]*x[i_2,])%*%beta.old), "+")
        x_p = x[i_1,]
        x_m = x[i_2,]

        derv.m <- -(x_p[rep(1:nrow(x_p), times=rep(length(i_2),nrow(x_p))),] - x_m[rep(1:length(i_2), length(i_1)),])
        grad.v <- apply( (exp(-c(t(m.ij)))/(1+exp(-c(t(m.ij))))) * derv.m, 2, mean) + pen.f(beta.c, a, lambda)*sign(beta.c)
        prox.el <- beta.old - grad.v*t.s
        G.t <- (beta.old-S.f(prox.el, lambda, t.s))/t.s
        beta.new <- beta.old - t.s*G.t

        if ( max(abs(beta.new - beta.old)) < eps ) break
        beta.old <- beta.new
      }

      if(max(abs(beta.new - beta.c)) < eps) break
      beta.c <- beta.new
    }
  }

  val <- structure(list(beta=round(beta.new, 4), lambda=lambda, penalty=penalty), class="auc.;pg.pen")

}

