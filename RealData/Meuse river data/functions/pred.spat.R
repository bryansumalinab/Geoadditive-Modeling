pred.spat <- function(model, coord0){
  K <- model$covs$K
  covs <- model$covs
  v_mode <- model$v_mode
  xi_estim <- model$xi_estim
  Sigma <- model$Sigma

  w1.0 <- coord0[,1]
  w2.0 <- coord0[,2]
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  
  xi_spat <- tail(xi_estim, K + 2)
  fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)
  
  # #####################
  ind.spat <- tail(1:length(xi_estim), K + 2)
  covarxi_spat <- Sigma[ind.spat, ind.spat]
  sd.spat <- sqrt(diag(Cspat.cent %*% covarxi_spat %*% t(Cspat.cent)))
  quantiles_spat <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, fit.w0, sd.spat)
  quantiles_w0 <- t(quantiles_spat)
  
  output <- data.frame("fit.spat" = fit.w0,
                       "lowerCI" = t(quantiles_spat)[,1],
                       "upperCI" = t(quantiles_spat)[,2])
  
  return(output)
}