Krig.Pois <- function(data, cov.method, knots = knots){
  
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  y <- data$y
  n <- length(y)
  xsmooth <- data$xsmooth
  coord <- cbind(data$w1, data$w2)
  xlin <- data$xlin
  X_linear <- cbind(1, xlin)
  
  B <- 30
  sm <- fsmooth.bspline(xsmooth, B)
  nbasis <- sm$nbasis
  Csmooth <- sm$Z
  
  # Spatial component
  knots <- knots
  K <- dim(knots)[1]
  covs <- cov.spatial(coord = coord,
                      K = K,
                      cov.method = cov.method,
                      knots = knots)
  
  C <- function(v) Matrix::Matrix(cbind(X_linear, Csmooth, covs$X,covs$Zw(v)), sparse = T)
  Q <- function(v) {
    result_matrix <- Matrix::bdiag(Matrix::Diagonal(n = dim(X_linear)[2], x = zeta), 
                                   exp(v[1]) * sm$G,
                                   Matrix::Diagonal(n = 2, x = zeta),
                                   exp(v[length(v) - 1]) * covs$Covw(v))
    return(result_matrix)
  }
  
  dimxi <- K + 2 + dim(X_linear)[2] + dim(Csmooth)[2]
  rk <- c(dim(sm$Z)[2], K)
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a <- b <- 1e-05
  # nu prior parameter for the penalty
  nu <- 3
  # Prior precision fixed effect
  zeta <- 1e-05 
  
  
  # Log conditional posterior of xi given v
  logpxi <- function(xi, Cv, Qv) {
    Cvxi <- as.numeric(Cv %*% xi)
    value <- sum(y * Cvxi - exp(Cvxi)) - .5 * t(xi) %*% Qv %*% xi
    as.numeric(value)
  }
  
  # Gradient of xi
  grad.logpxi <- function(xi,Cv, Qv){
    value <- Matrix::t(Cv)%*%(y - exp(as.numeric(Cv %*% xi))) - Qv%*%xi
    as.numeric(value)
  }
  
  Hess.logpxi <- function(xi,Cv, Qv){
    value <- - Matrix::t(Cv)%*%Matrix::Diagonal(x = exp(as.numeric(Cv %*% xi)))%*%Cv - Qv + Matrix::Diagonal(n = dimxi, x = 1e-08)
    value
  }
  
  # Laplace approximation to conditional posterior of xi
  Laplace <- function(xi0, Cv, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * solve(Hess.logpxi(xi0, Cv, Qv),
                                     grad.logpxi(xi0, Cv, Qv)))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- logpxi(xi0, Cv, Qv)
      while (logpxi(xi.new, Cv, Qv) <= logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    xistar <- xi0 
    return(xistar)
  }
  
  dist.ed <- function(w1, w2) {
    distances <- sqrt((w1 - w2)^2)
    return(distances)
  }
  dist <- dist.ed(coord[,1], coord[,2])

  v_init <- c(6, 0, log(1/max(covs$Z.dist)))
  Cv_init <- C(v_init)
  Qv_init <- Q(v_init)
  
  # Initial estimate for xi
  xi_hat <- Laplace(xi0 = rep(0,dimxi), 
                    Cv = Cv_init, 
                    Qv = Qv_init)
  log_pv <- function(v){
    vs <- v[-length(v)]
    vphi <- v[length(v)]
    Cv <- C(v)
    Cvxi <- as.numeric(Cv %*% xi_hat)
    Qv <- Q(v)
    e1 <- eigen(covs$Covw(v),only.values = T)$values
    e2 <- eigen(-Matrix::t(Cv)%*%Matrix::Diagonal(x = exp(Cvxi))%*%Cv + Qv,
                only.values = T)$values
    e1 <- e1[e1>0]
    e2 <- e2[e2>0]
    
    a1 <- sum(y * Cvxi - exp(Cvxi)) - .5 * t(xi_hat) %*% Qv %*% xi_hat
    a2 <- 0.5*sum(sapply(e1, log))
    a3 <- -0.5*sum(sapply(e2, log))
    a4 <- sum(0.5 * (rk + nu) * vs)
    a5 <- - sum(((0.5 * nu) + a) *  log(0.5*(nu * exp(vs)) + b))
    a6 <- (0.5 * nu) * vphi - ((0.5 * nu) + a) *  log(0.5*(nu * exp(vphi)) + b)
    value <- a1+a2+a3+a4+a5+a6
    return(as.numeric(value))
  }
  
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method = "Nelder-Mead", 
                  control = list(fnscale = -1, 
                                 reltol = 1e-5))$par
  
  Cv_mode <- C(v_mode)
  Qv_mode <- Q(v_mode)
  
  # Mode a posteriori estimate for xi
  xi_estim <- Laplace(xi_hat, 
                      Cv = Cv_mode, 
                      Qv = Qv_mode)
  
  Sigma <- -solve(Hess.logpxi(xi = xi_estim,
                              Cv = Cv_mode,
                              Qv = Qv_mode))
  output <- list("xi_estim" = xi_estim,
                 "Sigma" = Sigma,
                 "covs" = covs, 
                 "v_mode" = v_mode, 
                 "K" = K,
                 "sm" = sm) 
}