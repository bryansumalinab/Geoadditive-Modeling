rm(list = ls())

#Source the functions
file_list <- list.files(path = "functions", full.names = TRUE)
for (file in file_list) {
  source(file)
}

# Covariance functions:  circular, exponential, Matern, and spherical
cov.method <- "circular"

# Function s2
smooth_spat <- function(x, y) {
  return((x^3 + x*y + y^2)/25)
}
minysp <- -3
maxysp <- 3
n <- 1000 # sample size
sigma2 <- 0.10
beta0 <- 3

set.seed(0519)
knot.sp <- fields::cover.design(R = cbind(runif(n = n, min = minysp, max = maxysp), 
                                          runif(n = n, min = minysp, max = maxysp)), 
                                nd = 150)
knots <- as.matrix(knot.sp$design)
fsm <- function(x) cos(2*pi*x)
N <- 250 # number of simulations

CI_sm <- RBias_sm <- numeric(N)
CI_spat <- RBias_spat <- numeric(N)
PI_y0 <- numeric(N)
BIAS_mu <- RBias_mu <- numeric(N)

vmat <- matrix(NA, nrow = 3, ncol = N)
beta1 <- numeric(N)
varbeta1 <- numeric(N)
taumean <- numeric(N)

set.seed(0519)
for (ii in 1:N) {
  tryCatch({
    
    w1 <- runif(n, min = minysp, max = maxysp)
    w2 <- runif(n, min = minysp, max = maxysp)
    xlin <- runif(n)
    xsmooth <- runif(n, -1, 1)
    fsmooth <- cos(2 * pi * xsmooth)
    
    er <- rnorm(n, mean = 0, sd = sqrt(sigma2))
    ysp <- smooth_spat(x = w1, y = w2)
    y <- beta0 - 0.5 * xlin + fsmooth + ysp + er
    ymean <- beta0 - 0.5 * xlin + fsmooth + ysp
    
    dat <- data.frame(y, xlin, xsmooth, "w1" = w1, "w2" = w2)
    
    model <- Krig.Gauss(data = dat,
                        cov.method = cov.method,
                        knots = knots)
    
    nbasis <- model$sm$nbasis
    xi_estim <- model$xi_estim
    const.sm <- model$sm$const
    xsmooth <- model$sm$x
    B <- nbasis
    Sigma <- model$Sigma
    tau <- model$tau
    covs <- model$covs
    v_mode <- model$v_mode
    K <- model$covs$K
    
    fsm <- function(x) cos(2*pi*x)
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
    
    Sigma_smooth <- Sigma[ind, ind]
    covarxi_smooth <- mean(1/tau) * Sigma_smooth
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
    
    ################# spatial component ###########
    ################################################
    
    x <- seq(minysp, maxysp, length.out = 50)
    y <- seq(minysp , maxysp, length.out = 50)
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
    
    #####################
    ind.spat <- tail(1:length(xi_estim), K + 2)
    Sigma_spat <- Sigma[ind.spat, ind.spat]
    covarxi_spat <- mean(1/tau) * Sigma_spat
    sd.spat <- sqrt(diag(Cspat.cent %*% covarxi_spat %*% t(Cspat.cent)))
    quantiles_spat <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, fit.w0, sd.spat)
    quantiles_w0 <- t(quantiles_spat)
    
    #######
    ci_spat <- numeric(length(ysp.0))
    for (i in 1:length(ysp.0)) {
      ci_spat[i] <- as.numeric(ifelse(ysp.0[i] >= quantiles_w0[i, 1] & ysp.0[i] <= quantiles_w0[i, 2], 1, 0))
    }
    
    CI_spat[ii] <- mean(ci_spat)*100
    RBias_spat[ii] <- median(abs((ysp.0 - fit.w0)/ysp.0)*100)
    
    #############################################################
    ########### y prediction ###########
    
    ################# spatial component ###########
    ################################################
    
    x <- seq(minysp, maxysp, length.out = 10)
    y <- seq(minysp , maxysp, length.out = 10)
    grid <- expand.grid(x, y)
    w1.0 <- grid[,1]
    w2.0 <- grid[,2]
    
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
    
    xi_spat <- tail(xi_estim, K + 2)
    fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)
    
    ########################################################
    ########## Mean model prediction ###################
    ysp.0 <- smooth_spat(x = w1.0, y = w2.0)
    er.0 <- rnorm(n = ngrid, mean = 0, sd = sqrt(sigma2))
    xlin.0 <- runif(ngrid)
    mu.0 <- beta0 - 0.5*xlin.0 + fsm(xgrid) + ysp.0
    y0 <- mu.0 + er.0
    
    C0 <- cbind(1, xlin.0, Zsm.grid, Cspat.cent)
    fit.y0 <- as.numeric(C0 %*% xi_estim)
    
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
    
    #######
    ci_y0 <- numeric(length(y0))
    for (i in 1:length(y0)) {
      ci_y0[i] <- as.numeric(ifelse(y0[i] >= quantiles_y0[1, i] & y0[i] <= quantiles_y0[2, i], 1, 0))
    }
    
    ###### Prediction interval for y0
    PI_y0[ii] <- mean(ci_y0)*100
    
    ##### Bias and relative bias for mu
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
