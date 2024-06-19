plot.smooth <- function(model, term = NULL, text_size = 1.2, axis_size = 1.0) {
  
  if (!require("Matrix", character.only = TRUE)) {
  message(paste("Package Matrix", "is required but not installed."))
  }
  
  varnames <- colnames(model$X_smooth)
  sm.list <- model$sm.list
  nl <- dim(model$X_linear)[2]
  ns <- dim(model$X_smooth)[2]
  xi_estim <- model$xi_estim
  Sigma <- model$Sigma
  
  if (!is.null(term)) {
    if (term < 1 || term > ns) {
      stop("The specified smooth term index is out of range.")
    }
    plot_term <- term
  } else {
    stop("The specified smooth term index should be specified.")
  }
  smooth <- sm.list[[plot_term]]
  nbasis <- smooth$nbasis
  ngrid <- 500
  start <- nl + 1
  for (i in 1:plot_term){
  stop <- start + nbasis - 1
  ind <- start:stop
  start <- stop + 1
  }
  x0 <- seq(min(smooth$x), max(smooth$x), length.out = ngrid)
  const <- smooth$const
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

  plot(x0, fshat, type = "l", ylim = range(quantiles_y0),
       xlab = paste0(varnames[i]),
       ylab = paste0("f(", varnames[i], ")"), lwd = 2,
       cex.lab = text_size, cex.axis = axis_size)
  lines(x0, t(quantiles_y0)[,1], col = "red", lwd = 2)
  lines(x0, t(quantiles_y0)[,2], col = "red", lwd = 2)
  abline(h = 0, lty = 3, lwd = 1.5)
}
