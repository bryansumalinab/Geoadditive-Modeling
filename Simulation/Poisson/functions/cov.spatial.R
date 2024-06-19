cov.spatial <- function(coord, K, cov.method = "matern", knots){
  w1 <- coord[,1]
  w2 <- coord[,2]
  n <- length(w1)
  knots <- knots
  
  constX <- colMeans(cbind(w1, w2))
  X <- as.matrix(cbind(w1, w2) - matrix(rep(constX, each = n, byrow=T), nrow=n))
  
  k.dist <- matrix(0, K, K)
  k.dist[lower.tri(k.dist)] <- dist(knots)
  k.dist <- k.dist + t(k.dist)
  
  Z1.dist <- outer(w1,knots[,1],"-")
  Z2.dist <- outer(w2,knots[,2],"-")
  Z.dist <- sqrt(Z1.dist^2 + Z2.dist^2)
  
  if (cov.method == "exponential"){
    Scov <- function(x, v){
      vphi <- v[length(v)]
      x <- as.matrix(x)
      num.row <- nrow(x)
      num.col <- ncol(x)
      x <- as.vector(x)
      cov <- exp(-exp(vphi)*x)
      cov.mat<- matrix(cov,num.row,num.col)
      return(cov.mat)
    }
  }
  else if (cov.method == "gaussian")
    Scov <- function(x, v){
      vphi <- v[length(v)]
      x <- as.matrix(x)
      num.row <- nrow(x)
      num.col <- ncol(x)
      x <- as.vector(x)
      cov <- exp(-(exp(vphi)*x)^2)
      cov.mat<- matrix(cov,num.row,num.col)
      return(cov.mat)
    }
  else if (cov.method == "matern"){
    Scov <- function(x, v){
      vphi <- v[length(v)]
      x <- as.matrix(x)
      num.row <- nrow(x)
      num.col <- ncol(x)
      x <- as.vector(x)
      cov <- exp(-exp(vphi)*x)*(1 + exp(vphi)*x)
      cov.mat<- matrix(cov,num.row,num.col)
      return(cov.mat)
    }
  }
  else if (cov.method == "spherical"){
    Scov <- function(x, v){
      vphi <- v[length(v)]
      x <- as.matrix(x)
      num.row <- nrow(x)
      num.col <- ncol(x)
      x <- as.vector(x)
      cov <- ifelse(x <= (1/exp(vphi)),
                    1 - 1.5*exp(vphi)*x + 0.5*(exp(vphi)*x)^3,
                    0)
      cov.mat<- matrix(cov,num.row,num.col)
      return(cov.mat)
    }
  }
  
  else if (cov.method == "circular"){
    Scov <- function(x, v){
      vphi <- v[length(v)]
      x <- as.matrix(x)
      num.row <- nrow(x)
      num.col <- ncol(x)
      x <- as.vector(x)
      const <- pmin(x*exp(vphi), 1)
      cov <- 1 - (2/pi)*(const*sqrt(1 - const^2) + asin(const))
      cov.mat<- matrix(cov,num.row,num.col)
      return(cov.mat)
    }
  }
  constZ <- function(v) colMeans(Scov(x = Z.dist, v = v))
  Zw <- function(v) as.matrix(Scov(x = Z.dist, v = v) -
                                matrix(rep(constZ(v), each = n, byrow=T), nrow=n))
  
  Covw <- function(v) as.matrix(Scov(x = k.dist, v = v))
  output <- list("X" = X, "Zw" = Zw,
                 "Covw" = Covw, "Scov" = Scov,
                 "knots.s" = knots, 
                 "constX" = constX, 
                 "constZ" = constZ,
                 "K" = K,
                 "Z.dist" = Z.dist)
  
  return(output)
}