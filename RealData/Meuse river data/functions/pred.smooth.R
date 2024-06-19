pred.smooth <- function(model, ngrid = 500){
  output <- list(results = list())  # Initialize dataframes and plots lists
  
  varnames <- colnames(model$X_smooth)
  sm.list <- model$sm.list
  nl <- dim(model$X_linear)[2]
  ns <- dim(model$X_smooth)[2]
  xi_estim <- model$xi_estim
  Sigma <- model$Sigma
  
  start <- nl + 1
  for(i in 1:ns){
    smooth <- sm.list[[i]]
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
    sd <- sqrt(diag(Xsm %*% covarxi %*% t(Xsm)))
    quantiles_y0 <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, fshat, sd)
    
    output$results[[i]] <- data.frame(setNames(list(fshat,
                                                       t(quantiles_y0)[,1],
                                                       t(quantiles_y0)[,2]),
                                                  c(paste0("f(", varnames[i], ")"),
                                                    "lowerCI",
                                                    "upperCI")))

    start <- stop + 1
  }
  
  return(output)
}
