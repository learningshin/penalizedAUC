
pen.f <- function(x, a, lambda){
  val <- -lambda*(1-pmax(a*lambda-abs(x),0)/((a-1)*lambda)) *
    (abs(x)>lambda)
  return(val)
}

S.f <- function(x,lambda,t){
  val <- c(x>lambda*t)*(x-lambda*t) + c(abs(x)<=lambda*t)*0 + c(x< -lambda*t)*(x+lambda*t)
  return(val)
}
