Krig.Gauss <- function(data, cov.method, knots){
  
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  y <- data$y
  n <- length(y)
  coord <- cbind(data$w1, data$w2)
  knots <- knots
  K <- dim(knots)[1]
  covs <- cov.spatial(coord = coord,
                      K = K,
                      cov.method = cov.method,
                      knots = knots)
  
  C <- function(v) as.matrix(cbind(1,covs$X,covs$Zw(v)))
  Q <- function(v) {
    result_matrix <- Matrix::bdiag(Matrix::Diagonal(n = 3, x = zeta),
                                   exp(v[length(v) - 1]) * covs$Covw(v))
    return(result_matrix)
  }
  dimxi <- K + 3
  rk <- K
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a <- b <- 1e-05
  # nu prior parameter for the penalty
  nu <- 3
  # Prior precision fixed effect
  zeta <- 1e-05 
  
  # Mode a posteriori estimate of xi
  xihat <- function(v){
    Cv <- C(v)
    solve(Matrix::t(Cv)%*%Cv + Q(v) + Matrix::Diagonal(n = dimxi, x = 1e-10))%*%Matrix::t(Cv)%*%y
  }
  v_init <- c(1, log(1/max(covs$Z.dist)))
  xiv <- xihat(v_init)
  
  # Log - posterior of penalty vector v and and range parameter phi
  log_pv <- function(v){
    
    Cv <- C(v)
    Qv <- Q(v)
    CvC <- crossprod(Cv)
    vs <- v[-length(v)]
    vphi <- v[length(v)]
    
    e1 <- eigen(covs$Covw(v),only.values = T)$values
    e2 <- eigen(CvC + Qv,only.values = T)$values
    e1 <- e1[e1>0]
    e2 <- e2[e2>0]
    
    a1 <- -(0.5*n)*log(sum((y - Cv%*%xiv)^2) + t(xiv)%*% Qv %*% xiv)
    a2 <- 0.5*sum(sapply(e1, log))
    a3 <- -0.5*sum(sapply(e2, log))
    a4 <- sum(0.5 * (rk + nu) * vs)
    a5 <- - sum(((0.5 * nu) + a) *  log(0.5*(nu * exp(vs)) + b))
    a6 <- (0.5 * nu) * vphi - ((0.5 * nu) + a) *  log(0.5*(nu * exp(vphi)) + b)
    
    value <- a1 + a2 + a3 + a4 + a5 + a6
    
    return(as.numeric(value))
  }
  
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method = "Nelder-Mead", 
                  control = list(fnscale = -1, 
                                 reltol = 1e-3))$par
  
  phi <- 1/exp(v_mode[length(v_mode)])
  
  
  xi_estim <- xihat(v_mode)
  Cv <- C(v_mode)
  Qv <- Q(v_mode)
  
  Sigma <- solve(Matrix::t(Cv) %*% Cv + Qv + Matrix::Diagonal(n = dimxi, x = 1e-10))
  
  tau <- rgamma(n = 500, shape = n/2, 
                rate = as.numeric(0.5*(sum((y - Cv%*%xi_estim)^2)
                                       + t(xi_estim) %*% Qv %*% xi_estim)))
  
  output <- list("xi_estim" = xi_estim,
                 "Sigma" = Sigma,
                 "tau" = tau,
                 "covs" = covs, 
                 "v_mode" = v_mode, 
                 "K" = K)
}