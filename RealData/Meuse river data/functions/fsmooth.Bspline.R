fsmooth.bspline <- function(x, B){
  n <- length(x)
  Z <- cubicbs(x,lower = min(x),upper = max(x),K = B)$Bmatrix
  const <- colMeans(Z)
  Z <- Z - matrix(rep(const, each = n, byrow=T), nrow=n)
  
  penorder <- 2
  D <- Matrix::Diagonal(n = B)
  for (k in 1:penorder) D <- Matrix::diff(D)
  P <- Matrix::t(D) %*% D
  G <- P + Matrix::Diagonal(n = B, x = 1e-12)

  output <- list("Z" = Z,
                 "G" = G,
                 "x" = x,
                 "nbasis" = B,
                 "const" = const)
  
  return(output)
  
}
