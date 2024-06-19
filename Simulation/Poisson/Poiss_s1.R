rm(list = ls())

#Source the functions
file_list <- list.files(path = "functions", full.names = TRUE)
for (file in file_list) {
  source(file)
}

# Covariance functions:  circular, exponential, Matern, and spherical
cov.method <- "circular"

# Function s1
smooth_spat <- function(x, y) {
  z <- 0.5 - (x^2 + y^2)/(18)
  return(z)
}
minysp <- -3
maxysp <- 3
n <- 1000 # sample size
sderr <- 0.25
beta0 <- 3

set.seed(0519)
knot.sp <- fields::cover.design(R = cbind(runif(n = n, min = minysp, max = maxysp), 
                                          runif(n = n, min = minysp, max = maxysp)), 
                                nd = 150)
knots <- as.matrix(knot.sp$design)
fsm <- function(x) cos(2*pi*x)

N <- 250 # number of simulations

CI_sm <- RBias_sm <- numeric(N)
CI_spat <- RBias_spat <-  numeric(N)
PI_y0 <- numeric(N)
BIAS_mu <- RBias_mu <- numeric(N)

set.seed(0519)
for (ii in 1:N) {
  tryCatch({
    w1 <- runif(n = n, min = minysp, max = maxysp)
    w2 <- runif(n = n, min = minysp, max = maxysp)
    gridw <- cbind(w1, w2)
    
    xlin <- runif(n)
    xsmooth <- runif(n, -1, 1)
    fsmooth <- cos(2 * pi * xsmooth)
    
    ysp <- smooth_spat(x = w1, y = w2)
    er <- rnorm(n = n, mean = 0, sd = sderr)
    mu <- exp(beta0 - 0.5*xlin + fsmooth + ysp + er)
    y <- rpois(n = n, lambda = mu)
    coord <- gridw
    dat <- data.frame(y, mu, ysp, xlin, xsmooth, "w1" = coord[,1], "w2" = coord[,2], er)
    
    
    model <- Krig.Pois(data = dat,
                       cov.method = cov.method,
                       knots = knots)
    nbasis <- model$sm$nbasis
    xi_estim <- model$xi_estim
    const.sm <- model$sm$const
    xsmooth <- model$sm$x
    B <- nbasis
    Sigma <- model$Sigma
    covs <- model$covs
    v_mode <- model$v_mode
    K <- model$covs$K
    
    ################# (a) smooth component ###########
    ################################################
    
    ngrid <- 100
    xgrid <- seq(min(xsmooth), max(xsmooth), length.out = ngrid)
    ind <- (3): (2 + nbasis) 
    xi_smooth <- xi_estim[ind]
    
    Zsm.grid <- fsmooth.pred.bspline(x0 = xgrid,
                                     B = nbasis,
                                     const = const.sm,
                                     xmin = min(xsmooth),
                                     xmax = max(xsmooth))
    fshat.xgrid <- as.numeric(Zsm.grid %*% xi_smooth)
    
    covarxi_smooth <- Sigma[ind, ind]
    sd <- sqrt(diag(Zsm.grid %*% covarxi_smooth %*% t(Zsm.grid)))
    quantiles_xgrid <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, fshat.xgrid, sd)
    quantile_xsmooth <- t(quantiles_xgrid)
    fsxgrid <- fsm(xgrid) - mean(fsm(xgrid))
    
    
    ci_sm <- numeric(length(xgrid))
    for (i in 1:length(xgrid)) {
      ci_sm[i] <- as.numeric(ifelse(fsxgrid[i] >= quantiles_xgrid[1, i] & fsxgrid[i] <= quantiles_xgrid[2, i], 1, 0))
    }
    
    CI_sm[ii] <- mean(ci_sm)*100
    RBias_sm[ii] <- median(abs((fsxgrid - fshat.xgrid)/fsxgrid)*100)
    
    
    ################# (b) spatial component ###########
    ################################################
    
    x <- seq(minysp, maxysp, length.out = 50)
    y <- seq(minysp, maxysp, length.out = 50)
    grid <- expand.grid(x, y)
    w1.0 <- grid[,1]
    w2.0 <- grid[,2]
    z <- smooth_spat(x = w1.0, y = w2.0)
    ysp.0 <- z - mean(z)
    
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
    
    xi_spat <- tail(xi_estim, K + 2)
    fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)
    
    ind.spat <- tail(1:length(xi_estim), K + 2)
    covarxi_spat <- Sigma[ind.spat, ind.spat]
    sd.spat <- sqrt(diag(Cspat.cent %*% covarxi_spat %*% t(Cspat.cent)))
    quantiles_spat <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, fit.w0, sd.spat)
    quantiles_w0 <- t(quantiles_spat)
    
    ci_spat <- numeric(length(ysp.0))
    for (i in 1:length(ysp.0)) {
      ci_spat[i] <- as.numeric(ifelse(ysp.0[i] >= quantiles_w0[i, 1] & ysp.0[i] <= quantiles_w0[i, 2], 1, 0))
    }
    CI_spat[ii] <- mean(ci_spat)*100
    RBias_spat[ii] <- median(abs((ysp.0 - fit.w0)/ysp.0)*100)
    
    #############################################################
    ########### y prediction ###########
    
    ###### spatial component
    
    x <- seq(minysp, maxysp, length.out = 10)
    y <- seq(minysp, maxysp, length.out = 10)
    grid <- expand.grid(x, y)
    w1.0 <- grid[,1]
    w2.0 <- grid[,2]
    z <- smooth_spat(x = w1.0, y = w2.0)
    
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
    
    xi_spat <- tail(xi_estim, K + 2)
    fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)
    
    ########################################################
    er.0 <- rnorm(n = ngrid, mean = 0, sd = sderr)
    xlin.0 <- runif(ngrid)
    fsm.0 <- fsm(xgrid)
    mu.0 <- exp(beta0 - 0.5*xlin.0 + fsm.0 + z)
    y0 <- rpois(n = ngrid, 
                lambda = exp(beta0 - 0.5*xlin.0
                             + fsm.0 + z + er.0))
    
    C0 <- cbind(1, xlin.0, Zsm.grid, Cspat.cent)
    fit.y0 <- exp(C0 %*% xi_estim)
    
    # Compute mean and variance for log(mu)
    logmu <- C0%*%xi_estim
    logmu.var <- diag(C0%*%Sigma%*%t(C0))
    
    
    # Generate Poisson samples
    r.pois <- list()
    N <- 1000
    for (i in 1:length(logmu)) {
      rn <- rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
      mu <- exp(rn)
      r.pois[[i]] <- rpois(N,lambda = mu)
    }
    quantiles_y0 <- sapply(r.pois, function(x) quantile(x, probs = c(0.025, 0.975)))
    
    ci_y0 <- numeric(length(y0))
    for (i in 1:length(y0)) {
      ci_y0[i] <- as.numeric(ifelse(y0[i] >= quantiles_y0[1, i] & y0[i] <= quantiles_y0[2, i], 1, 0))
    }
    
    ###### Prediction interval for y0
    PI_y0[ii] <- mean(ci_y0)*100
    
    ##### Bias and reltive bias for mu
    BIAS_mu[ii] <- mean(mu.0 - fit.y0)
    RBias_mu[ii] <- median(abs((mu.0 - fit.y0)/mu.0)*100)
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(ii)
}


# results
result <- data.frame(
  Metric = c("Bias", "Relative Bias", "Coverage"),
  smooth = c(NA, mean(RBias_sm), mean(CI_sm)),
  spatial = c(NA, mean(RBias_spat), mean(CI_spat)),
  mu = c(mean(BIAS_mu), mean(RBias_mu), NA),
  y0 = c(NA, NA, mean(PI_y0))
)

result
