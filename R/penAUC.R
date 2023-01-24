
penAUC <- function(x, y, loss=c("exponential", "logistic"), penalty=c("lasso","ridge","scad"),
                   lambda, learn.rate=0.01, eps=1e-3, max.iter=1e3, a=3.7){

  # Response +1 or -1
  if( sum(unique(y) %in% c(-1,1)) != 2 )
    warning(print(" The response y shoud take the value +1 or -1."))

  # Specify the loss function and penalty function
  if( length(loss) != 1 )
    warning(print(" Loss function should be assigned as either 'exponential' or 'logistic'."))

  if(  length(penalty) !=1 )
    warning(print(" Penalty function should be applied, 'Lasso', 'Ridge' or 'SCAD'."))

  if( is.null(lambda) )
    warning(print(" The regularization parameter should be specified"))

  loss <- match.arg(loss)
  penalty <- match.arg(penalty)

  # Fit penalizedAUC
  if(loss == "exponential"){
    obj <- auc.exp.pen(x, y, penalty, lambda=lambda, eps=eps, max.iter=max.iter, learn.rate=learn.rate, a=a)$beta
  } else{
    obj <- auc.log.pen(x, y, penalty, lambda=lambda, eps=eps, max.iter=max.iter, learn.rate=learn.rate, a=a)$beta
  }
  val <- structure(list(beta=round(matrix(obj,ncol=1),4), lambda=lambda, penalty = penalty, loss=loss), class="penAUC")
}

