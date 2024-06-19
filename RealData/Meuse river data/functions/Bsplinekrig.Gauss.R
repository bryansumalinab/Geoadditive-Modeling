Bsplinekrig.Gauss <-  function(y, X_linear = NULL, 
                                   X_smooth = NULL, coord,
                                   knots, cov.method){
  
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  y <- y
  n <- length(y)
  coord <- coord
  cov.method <- cov.method
  
  if (!is.null(X_linear)) {
    X_linear <- model.matrix(~ ., data = X_linear)
  } else {
    X_linear <- matrix(1, nrow = length(y))
    colnames(X_linear) <- "(Intercept)"
  }
  nl <- ncol(X_linear)
  
  
  if (!is.null(X_smooth)){
    B <- 30 # Number of B-spline basis
    ns <- ncol(X_smooth) # Number of smooth terms
    
    sm.list <- list()
    for (i in 1:ns) {
      x <- X_smooth[,i]
      sm.list[[i]] <- fsmooth.bspline(x, B = B)
    }
    
    Csmooth <- Matrix::Matrix(do.call(cbind, lapply(sm.list, function(sm) sm$Z)), sparse = TRUE)
    Qsmooth <- function(v) Matrix::bdiag(
      lapply(1:ns, function(i) {
        exp(v[i]) * sm.list[[i]]$G
      }))
    
    # spatial component
    if (missing(knots)) {
      if (!require("fields", character.only = TRUE)) {
        message(paste("Package fields", "is required but not installed."))
      }
      K <- max(20, min(length(w1)/4, 150))
      knot.sp <- fields::cover.design(R = cbind(coord[,1],coord[,2]), nd = K)
      knots <- as.matrix(knot.sp$design)
    }
    knots <- knots
    K <- dim(knots)[1]
    
    covs <- cov.spatial(coord = coord,
                        K = K,
                        cov.method = cov.method,
                        knots = knots)
    C <- function(v) Matrix::Matrix(cbind(X_linear, Csmooth, covs$X,covs$Zw(v)), sparse = T)
    
    Q <- function(v) {
      result_matrix <- Matrix::bdiag(Matrix::Diagonal(n = nl, x = zeta), 
                                     Qsmooth(v),
                                     Matrix::Diagonal(n = 2, x = zeta),
                                     exp(v[length(v) - 1]) * covs$Covw(v))
      return(result_matrix)
    }
    dimxi <- K + 2 + nl + dim(Csmooth)[2]
    rk <- c(sapply(sm.list, function(sm) dim(sm$Z)[2]), K)
  } else {
    ns <- 0
    B <- NULL
    X_smooth <- NULL
    sm.list <- NULL
    # spatial component
    if (missing(knots)) {
      K <- max(20, min(length(w1)/4, 150))
      knot.sp <- fields::cover.design(R = cbind(coord[,1],coord[,2]), nd = K)
      knots <- as.matrix(knot.sp$design)
    }
    knots <- knots
    K <- dim(knots)[1]
    
    covs <- cov.spatial(coord = coord,
                        K = K,
                        cov.method = cov.method,
                        knots = knots)
    C <- function(v) Matrix::Matrix(cbind(X_linear, covs$X,covs$Zw(v)), sparse = T)
    
    Q <- function(v) {
      result_matrix <- Matrix::bdiag(Matrix::Diagonal(n = nl + 2, x = zeta),
                                     exp(v[length(v) - 1]) * covs$Covw(v))
      return(result_matrix)
    }
    dimxi <- K + nl + 2
    rk <- K 
  }
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a <- b <- 1e-05
  # nu prior parameter for the penalty
  nu <- 3
  # Prior precision fixed effect
  zeta <- 1e-05 
  
  # Mean psoterior estimate of xi
  xihat <- function(v){
    Cv <- C(v)
    solve(Matrix::t(Cv)%*%Cv + Q(v) + Matrix::Diagonal(n = dimxi, x = 1e-10))%*%Matrix::t(Cv)%*%y
  }
  if(ns > 0){
    v_init <- c(rep(5, ns), 0, log(1/max(covs$Z.dist))) 
  } else {
    v_init <- c(0, log(1/max(covs$Z.dist))) 
  }
  xiv <- xihat(v_init)
  # Log - posterior of penalty vector v and and range parameter phi
  log_pv <- function(v){
    
    Cv <- C(v)
    Qv <- Q(v)
    CvC <- Matrix::crossprod(Cv)
    vs <- v[-length(v)]
    vphi <- v[length(v)]
    
    e1 <- eigen(covs$Covw(v),only.values = T)$values
    e2 <- eigen(CvC + Qv,only.values = T)$values
    e1 <- e1[e1>0]
    e2 <- e2[e2>0]
    
    a1 <- -(0.5*n)*log(sum((y - Cv%*%xiv)^2) + Matrix::t(xiv)%*% Qv %*% xiv)
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
  Sigma1 <- solve(Matrix::t(Cv) %*% Cv + Qv + Matrix::Diagonal(n = dimxi, x = 1e-10))
  
  tau <- rgamma(n = 500, shape = n/2, 
                rate = as.numeric(0.5*(sum((y - Cv%*%xi_estim)^2)
                                       + Matrix::t(xi_estim) %*% Qv %*% xi_estim)))
  Sigma <- mean(1/tau)*Sigma1
  
  fitted.values <- as.numeric(Cv%*%xi_estim)
  
  # Hessian of log-likelihood
  Hess.loglik <- mean(tau)*(t(Cv)%*%Cv)
  log_lik <- (0.5*n)*log(mean(tau)) - 0.5*mean(tau)*sum((y - Cv%*%xiv)^2)
  
  ED <- diag(Sigma %*% Hess.loglik)
  BIC <- as.numeric(-2*log_lik + sum(ED)*log(length(y)))
  
  
  ##### Test for smooth terms #############
  
  if(ns > 0){
    
    # Moore-Penrose inverse
    MPinv <- function(A) {
      svd_decomp <- svd(A)
      U <- svd_decomp$u
      D <- svd_decomp$d
      V <- svd_decomp$v
      D_inv <- diag(1/D, nrow = length(D), ncol = length(D))
      D_inv[is.infinite(D_inv)] <- 0
      result <- V %*% D_inv %*% t(U)
      return(result)
    }
    
    ED_smooth_result <- numeric(ns)
    Tstat_result <- numeric(ns)
    pval_result <- numeric(ns)
    
    start <- nl + 1
    
    ngrid <- 500
    for(ss in 1:ns){
      smooth <- sm.list[[ss]]
      nbasis <- smooth$nbasis
      x0 <- seq(min(smooth$x), max(smooth$x), length.out = ngrid)
      const <- smooth$const
      stop <- start + nbasis - 1
      ind <- start:stop
      
      xi_smooth <- xi_estim[ind]
      Xsm <- fsmooth.pred.bspline(x0 = x0,
                                  B = nbasis,
                                  const = const,
                                  xmin = min(smooth$x),
                                  xmax = max(smooth$x))
      fshat <- as.numeric(Xsm %*% xi_smooth)
      
      xi_estim_smooth <- xi_estim[ind]
      covarxi <- Sigma[ind, ind]
      V <- Xsm %*% covarxi %*% t(Xsm)
      Vinv <- MPinv(V)
      Tstat <- as.numeric(t(fshat) %*% Vinv %*% matrix(fshat))
      EDsmooth <- sum(ED[ind])
      pval <- stats::pgamma(Tstat, shape = EDsmooth/2, scale = 2, lower.tail = FALSE)
      ED_smooth_result[ss] <- EDsmooth
      Tstat_result[ss] <- Tstat
      pval_result[ss] <- pval
      
      start <- stop + 1
    }
    
    results_smooth <- data.frame(
      "edf" = ED_smooth_result,
      "Tr" = Tstat_result,
      "p value" = pval_result
    )
    rownames(results_smooth) <- colnames(X_smooth)
  } else {
    results_smooth <- NULL
  }
  
  ##### Linear terms #############
  
  beta_est <- xi_estim[1:nl]
  beta_sd <- sqrt(diag(Sigma)[1:nl])
  beta_CI <- matrix(0, nrow = nl, ncol = 2)
  for (ll in 1:nl) {
    beta_CI[ll, 1] <- beta_est[ll] - stats::qnorm(1- (1-0.95)/2) *
      beta_sd[ll]
    beta_CI[ll, 2] <- beta_est[ll] + stats::qnorm(1- (1-0.95)/2) *
      beta_sd[ll]
  }
  
  results_linear <- data.frame(
    "Estimate" = beta_est,
    "sd" = beta_sd,
    "CI95lower" = beta_CI[, 1],
    "CI95upper" = beta_CI[, 2]
  )
  rownames(results_linear) <- colnames(X_linear)
  
  output <- list("xi_estim" = xi_estim,
                 "Sigma1" = Sigma1,
                 "Sigma" = Sigma,
                 "covs" = covs,
                 "nbasis" = B, 
                 "v_mode" = v_mode, 
                 "K" = K,
                 "sm.list" = sm.list,
                 "X_linear" = X_linear,
                 "X_smooth" = X_smooth,
                 "coord" = coord,
                 "Cv_mode" = Cv,
                 "tau" = tau,
                 "fitted.values" = fitted.values,
                 "ED" = ED,
                 "BIC" = BIC,
                 "result_smooth" = results_smooth,
                 "results_linear" = results_linear)
}