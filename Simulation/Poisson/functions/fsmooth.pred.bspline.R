fsmooth.pred.bspline <- function(x0, xmin, xmax, B, const){
  Z <- cubicbs(x0,lower = xmin,upper = xmax,K = B)$Bmatrix
  Z <- Z - matrix(rep(const, each = length(x0), byrow=T), nrow=length(x0))
  return(Z)
}
