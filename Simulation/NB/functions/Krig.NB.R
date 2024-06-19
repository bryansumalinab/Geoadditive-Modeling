Krig.NB <- function(data, cov.method, knots){
  
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
                                   exp(v[length(v) - 2]) * covs$Covw(v))
    return(result_matrix)
  }
  
  dimxi <- K + 2 + dim(X_linear)[2] + dim(Csmooth)[2]
  rk <- c(dim(sm$Z)[2], K)
  
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a.delta <- 1e-05
  b.delta <- 1e-05
  # Prior for overdispersion parameter phi
  a.disp <- 1e-05
  b.disp <- 1e-05
  
  # nu prior parameter for the penalty
  nu <- 3
  # Prior precision fixed effect
  zeta <- 1e-05 
  
  
  # Negative binomial GLM with log-link
  mu.nb <- function(xi, v) exp(as.numeric(C(v) %*% xi))
  var.nb <- function(xi, v) {
    muval <- mu.nb(xi, v)
    res <- muval + (1 / exp(v[length(v) - 1])) * (muval ^ 2)
    return(res)
  }
  W.nb <- function(xi, v) {
    muval <- exp(as.numeric(C(v) %*% xi))
    varval <- muval + (1 / exp(v[length(v) - 1])) * (muval ^ 2)
    res <- Matrix::Diagonal(x = ((muval) ^ 2) * (1 / varval))
    return(res)
  }

  D.nb <- function(xi, v) Matrix::Diagonal(x = 1/mu.nb(xi, v))
  M.nb <- function(xi, v) Matrix::Diagonal(x = y - mu.nb(xi, v))
  V.nb <- function(xi, v){
    muval <- mu.nb(xi, v)
    varval <-  muval + (1 / exp(v[length(v) - 1])) * (muval ^ 2)
    res <- Matrix::Diagonal(x = muval * (1/varval - (muval/(varval^2)) * 
                                           (1 + 2*muval*(1/exp(v[length(v) - 1])))))
    return(res)
  }
  gamma.nb <- function(xi, v) {
    muval <- mu.nb(xi, v)
    res <- exp(v[length(v) - 1]) * log(muval / (muval + exp(v[length(v) - 1])))
    return(res)
  }
  bgamma.nb <- function(xi, v) - (exp(v[length(v) - 1])^2) *
    log(exp(v[length(v) - 1])/(exp(v[length(v) - 1]) + mu.nb(xi, v)))
  
  # Log conditional posterior of xi given v
  log_pxi <- function(xi, v) {
    value <- (1/exp(v[length(v) - 1])) * sum((y * gamma.nb(xi, v)) - bgamma.nb(xi, v)) -
      .5 * Matrix::t(xi) %*% Q(v) %*% xi
    return(as.numeric(value))
  }
  
  Grad.logpxi <- function(xi,v){
    Cv <- C(v)
    muval <- exp(as.numeric(Cv %*% xi))
    varval <- muval + (1 / exp(v[length(v) - 1])) * (muval ^ 2)
    W.nbval <- Matrix::Diagonal(x = ((muval) ^ 2) * (1 / varval))
    D.nbval <- Matrix::Diagonal(x = 1/muval)
    value <- as.numeric(Matrix::t(Cv)%*%W.nbval%*%D.nbval%*%(y - muval) -
                          Q(v)%*%xi)
    return(value)
  }
  
  # Hessian of parameter xi
  Hess.logpxi <- function(xi,v){
    Cv <- C(v)
    value <- Matrix::t(Cv)%*%M.nb(xi, v)%*%V.nb(xi, v)%*%Cv -
      Matrix::t(Cv)%*%W.nb(xi, v)%*%Cv - Q(v) + Matrix::Diagonal(n = dimxi, x = 1e-08)
    value
  }
  
  # Laplace approximation to conditional posterior of xi
  # using Newton-Raphson algorithm
  NR_xi <- function(xi0, v, Hess.logpxi, Grad.logpxi, log_pxi){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * Matrix::solve(Hess.logpxi(xi0, v),
                                     Grad.logpxi(xi0, v)))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- log_pxi(xi0, v)
      while (log_pxi(xi.new, v) <= logpxi.current) {
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

  v_init <- c(6, 0, 0, log(1/max(covs$Z.dist)))
  
  # Initial estimate for xi
  xi_hat <- NR_xi(xi0 = rep(0,dimxi), v = v_init, 
                  Hess.logpxi = Hess.logpxi, 
                  Grad.logpxi = Grad.logpxi, 
                  log_pxi = log_pxi)
  
  log_pv <- function(v){
    vs <- v[1:(length(v) - 2)]
    vphi <- v[length(v)]
    vdisp <- v[length(v) - 1]
    Cv <- C(v)
    Qv <- Q(v)
    Cvxi <- as.numeric(Cv %*% xi_hat)
    
    muval <- mu.nb(xi_hat, v)
    varval <-  muval + (1 / exp(vdisp)) * (muval ^ 2)
    Vnb <- Matrix::Diagonal(x = muval * (1/varval - (muval/(varval^2)) *
                                           (1 + 2*muval*(1/exp(vdisp)))))
    Wnb <- Matrix::Diagonal(x = ((muval) ^ 2) * (1 / varval))
    gammanb <- exp(vdisp) * log(muval / (muval + exp(vdisp)))
    bgammanb <- (-1) * (exp(vdisp)^2) * log(exp(vdisp)/(exp(vdisp) + muval))
    
    e1 <- eigen(covs$Covw(v),only.values = T)$values
    e2 <- eigen(-Matrix::t(Cv)%*%(M.nb(xi = xi_hat, v = v)%*%Vnb-Wnb)%*%Cv + Qv,
                only.values = T)$values
    e1 <- e1[e1>0]
    e2 <- e2[e2>0]
    
    a1 <- sum((1/exp(vdisp))*((y * gammanb) - bgammanb) +
                lgamma(y + exp(vdisp)) - lgamma(exp(vdisp))) - 0.5 * sum((xi_hat * Qv %*% xi_hat)) 
    a2 <- 0.5*sum(sapply(e1, log))
    a3 <- -0.5*sum(sapply(e2, log))
    a4 <- sum(0.5 * (rk + nu) * vs)
    a5 <- - sum(((0.5 * nu) + a.delta) *  log(0.5*(nu * exp(vs)) + b.delta))
    a6 <- (0.5 * nu) * vphi - ((0.5 * nu) + a.delta) *  log(0.5*(nu * exp(vphi)) + b.delta)
    a7 <- a.disp*vdisp - b.disp*exp(vdisp)
    value <- a1+a2+a3+a4+a5+a6+a7
    return(as.numeric(value))
  }
  
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method = "Nelder-Mead", 
                  control = list(fnscale = -1, 
                                 reltol = 1e-05))$par
  
  Cv_mode <- C(v_mode)
  Qv_mode <- Q(v_mode)
  
  # Mode a posteriori estimate for xi
  xi_estim <- NR_xi(xi0 = xi_hat, v = v_mode, 
                    Hess.logpxi = Hess.logpxi, 
                    Grad.logpxi = Grad.logpxi, 
                    log_pxi = log_pxi)
  
  Sigma <- -solve(Hess.logpxi(xi = xi_estim,
                              v = v_mode))
  
  output <- list("xi_estim" = xi_estim,
                 "Sigma" = Sigma,
                 "covs" = covs, 
                 "v_mode" = v_mode, 
                 "K" = K,
                 "sm" = sm) 
}