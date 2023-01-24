
cv.penAUC <- function(x, y, nfold=5, loss=c("exponential", "logistic"), penalty=c("lasso","ridge","scad"),
                      criteria=NULL, lambda.set=NULL, learn.rate=0.01, eps=1e-3, max.iter=1e3, a=3.7, seed=10){

  # Response +1 or -1
  if( sum(unique(y) %in% c(-1,1)) != 2 )
    warning(print(" The response y shoud take the value +1 or -1."))

  # Specify the loss function and penalty function
  if( length(loss) != 1 )
    warning(print(" Loss function should be assigned as either 'exponential' or 'logistic'."))

  if(  length(penalty) !=1 )
    warning(print(" Penalty function should be applied, 'lasso', 'ridge' or 'scad'."))

  if( is.null(lambda.set) ) lambda.set = 1.6^(-9:4)

  if( is.null(criteria) ){
    if(penalty != "scad"){
      criteria = "CV"
    } else{
      criteria = "AUCIC"
    }
  }

  loss <- match.arg(loss)
  penalty <- match.arg(penalty)
  len.lambda <- length(lambda.set)

  if( criteria == "CV" ){

    ### k-fold Cross validation AUC
    set.seed(seed)
    folds  <- sample(cut(seq(1,n),breaks=nfold,labels=FALSE),n)
    cv.result <- matrix(NA, nfold, len.lambda)
    p <- ncol(x)

    for(i in 1:nfold){

      ind = which(folds ==i)
      x.tr = x[-ind,]
      y.tr = y[-ind]
      x.ts = x[ind,]
      y.ts = y[ind]

      tr.i_p <- which(y.tr==1)
      tr.i_n <- which(y.tr==-1)

      ts.i_p <- which(y.ts==1)
      ts.i_n <- which(y.ts==-1)

      cv.obj <- matrix(NA,p,len.lambda)

      for(k in 1:len.lambda){
        lamb <- lambda.set[k]

        obj.tr <- penAUC(x.tr, y.tr, loss, penalty, lambda=lamb, learn.rate, eps, max.iter, a)
        cv.obj[,k]  <- obj.tr$beta
      }

      # test set for AUC
      fhat.ts <- as.matrix(x.ts) %*% cv.obj

      nm <- c(length(ts.i_p)*length(ts.i_n))
      if (length(ts.i_p)==0 || length(ts.i_n)==0)
        auc <- NA
      if (  (length(ts.i_n)!=1 && length(ts.i_p)!=1) ){
        auc <- rowSums(sapply(ts.i_p, function(v) apply(fhat.ts[v,] > t(fhat.ts[ts.i_n,]), 1, sum))) / nm
      } else{
        auc <- sum(sapply(ts.i_p, function(v) apply(fhat.ts[v] > t(fhat.ts[ts.i_n]), 1, sum))) / nm
      }
      cv.result[i,] <- auc
    }
    aver.auc <- apply(na.omit(cv.result), 2, mean)

    # choose lambda that maximum auc
    lamb.min.ind <- which.max(aver.auc)[1]
    result <- list(opt.lambda = round(lambda.set[lamb.min.ind],4), penalty = penalty, loss = loss,
                   criteria=criteria)

  } else if( (criteria == "AUCIC") & (penalty=="scad") ){

    obj.lambda.set = sapply(1:len.lambda, function(t) penAUC(x, y, loss, penalty, lambda=lambda.set[t], learn.rate, eps, max.iter, a)$beta)
    obj.lambda.set

    non.zero.scad = apply(obj.lambda.set!=0, 2, sum)
    cand.aucic = rep(NA,len.lambda)

    i_p <- which(y==1)
    i_m <- which(y==-1)

    for(m in 1:len.lambda){
      m.ij  <- outer( c((y[i_p]*x[i_p,])%*%obj.lambda.set[,m]), c((y[i_m]*x[i_m,])%*%obj.lambda.set[,m]), "+")
      cand.aucic[m] <- mean(exp(-m.ij)) + non.zero.scad[m]*log(n)/n
    }

    result <- list( opt.lambda = round(lambda.set[which.min(cand.aucic)],4), penalty = penalty, loss=loss,
                    criteria = criteria)
  }

  val <- structure(result, class="cv.penAUC")
}




# cv.obj = cv.penAUC(x,y,loss="exponential",penalty="scad",seed=6)
# cv.obj



