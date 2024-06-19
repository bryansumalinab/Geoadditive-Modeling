pred.Gauss <- function(model, X_linear0, 
                      X_smooth0, coord0,
                      knots, cov.method){
  
  n <- nrow(coord0)
  K <- model$covs$K
  covs <- model$covs
  v_mode <- model$v_mode
  xi_estim <- model$xi_estim
  Sigma <- model$Sigma1
  tau <- model$tau
  
  if (dim(model$X_linear)[2] > 1){
    if (missing(X_linear0)){
      stop("Missing linear predictor(s).")
    }
    if ((dim(model$X_linear[, -1, drop = FALSE])[2] !=  dim(X_linear0)[2])) {
      stop("The number of linear predictors are not the same in the original model.")
    }
    X_linear0 <- model.matrix(~ ., data = X_linear0)
  } else {
    X_linear0 <- matrix(1, nrow = n)
    colnames(X_linear0) <- "(Intercept)"
  }
  
  
  if (!is.null(model$X_smooth)){
    if (missing(X_smooth0)){
      stop("Missing smooth predictor(s).")
    }
    if (dim(model$X_smooth)[2] !=  dim(X_smooth0)[2]) {
      stop("The number of linear predictors are not the same in the original model.")
    }
    ns <- ncol(X_smooth0) # Number of smooth terms
    sm.list <- list()
    for (i in 1:ns) {
      x0 <- X_smooth0[,i]
      smooth <- model$sm.list[[i]]
      nbasis <- smooth$nbasis
      const <- smooth$const
      sm.list[[i]] <- fsmooth.pred.bspline(x0 = x0,
                                           B = nbasis,
                                           const = const,
                                           xmin = min(smooth$x),
                                           xmax = max(smooth$x))
    }
    Csmooth0 <- Matrix::Matrix(do.call(cbind, lapply(sm.list, function(sm) sm)), sparse = TRUE)
    
    w1.0 <- coord0[,1]
    w2.0 <- coord0[,2]
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat0 <- as.matrix(cbind(w1.0, w2.0, Zw.0) - 
                          matrix(rep(constspat, 
                                     each = length(w1.0), byrow=T), 
                                 nrow=length(w1.0)))
    C0 <- cbind(X_linear0, Csmooth0, Cspat0)
    fit.y0 <- as.numeric(C0 %*% xi_estim)
  } else {
    w1.0 <- coord0[,1]
    w2.0 <- coord0[,2]
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat0 <- as.matrix(cbind(w1.0, w2.0, Zw.0) - 
                          matrix(rep(constspat, 
                                     each = length(w1.0), byrow=T), 
                                 nrow=length(w1.0)))
    C0 <- cbind(X_linear0, Cspat0)
    fit.y0 <- as.numeric(C0 %*% xi_estim)
  }
  
  
  # Function to calculate mean prediction interval
  y_sim <- function(C0, tauvec, mu, D) {
    nsim <- length(tauvec)
    # Generate multivariate normal samples
    samples <- lapply(1:nsim, function(i) {
      covariance <- (1/tauvec[i]) * D
      thetas <- MASS::mvrnorm(n = 1, mu = mu, Sigma = covariance)
      n0 <- dim(C0)[1]
      
      values <- sapply(1:n0, function(j) {
        rnorm(n = 1, mean = as.numeric(C0[j,] %*% thetas), sd = sqrt(1/tauvec[i]))
      })
      return(values)
    })
    output <- do.call(rbind, samples)
    return(output)
  }
  
  ########### Posterior predictive samples
  y_samples <- y_sim(C0 = C0, mu = xi_estim,
                     tauvec = tau, D = Sigma)
  # Calculate quantiles for each column
  quantiles_y0 <- apply(y_samples, 2, quantile, c(0.025, 0.975))
  
  output <- data.frame("fit.y" = fit.y0,
                       "lowerCI" = quantiles_y0[1,],
                       "upperCI" = quantiles_y0[2,])
  
  return(output)
}
