rm(list = ls())

#Source the functions
file_list <- list.files(path = "functions", full.names = TRUE)
for (file in file_list) {
  source(file)
}

# Covariance functions:  circular, exponential, Matern, and spherical
cov.method <- "matern"

sill <- 0.5 # sill parameter
phi <- 0.15 # range parameter
sigma2 <- 0.10 # error variance
n <- 1000 # sample size

set.seed(0519)
# Two-dimensional knots
knot.sp <- fields::cover.design(R = cbind(runif(n = n, min = 0, max = 1),
                                          runif(n = n, min = 0, max = 1)),
                                nd = 150)
knots <- as.matrix(knot.sp$design)

N <- 250 # number of simulations
BIAS_mu <- RBias_mu <- rep(NA, N)
PI_y0 <- rep(NA, N)

set.seed(0519)
for (ii in 1:N) {
  tryCatch({
    
    w1 <- runif(n = n, min = 0, max = 1)
    w2 <- runif(n = n, min = 0, max = 1)
    gridw <- cbind(w1, w2)
    simdata <- geoR::grf(cov.pars = c(sill, phi),
                   grid = gridw, 
                   cov.model = cov.method, 
                   kappa = 1.5)
    
    er <- rnorm(n, mean = 0, sd = sqrt(sigma2))
    ysp <- simdata$data
    y <- 3.0 + ysp + er
    ymean <- 3.0 + ysp
    coord <- simdata$coords
    
    data <- data.frame(y, ymean,
                       "w1" = coord[,1],
                       "w2" = coord[,2])
    size.x0 <- floor(0.1 * nrow(data))
    ind.x0 <- sample(1:nrow(data), size = size.x0, replace = FALSE)
    # use for prediction
    data.pred <- data[ind.x0, ]
    # use for estimation
    data.est <- data[-ind.x0, ]
    
    model <- Krig.Gauss(data = data.est,
                        cov.method = cov.method,
                        knots = knots)
    
    xi_estim <- model$xi_estim
    Sigma <- model$Sigma
    tau <- model$tau
    covs <- model$covs
    v_mode <- model$v_mode
    K <- model$covs$K
    #############################################################
    ########### Prediction ###########
    w1.0 <- data.pred$w1
    w2.0 <- data.pred$w2
    y0 <- data.pred$y
    mu.0 <- data.pred$ymean
    
    Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
    Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
    Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
    
    Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
    constspat <- c(covs$constX, covs$constZ(v_mode))
    Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - 
                              matrix(rep(constspat, each = length(w1.0), byrow=T), 
                                     nrow=length(w1.0)))
    xi_spat <- tail(xi_estim, K + 2)
    C0 <- cbind(1, Cspat.cent)
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
    pi_y0 <- rep(NA, length(y0))
    for (i in 1:length(y0)) {
      pi_y0[i] <- as.numeric(ifelse(y0[i] >= quantiles_y0[1, i] & y0[i] <= quantiles_y0[2, i], 1, 0))
    }
    
    ###### Prediction interval for y0
    PI_y0[ii] <- mean(pi_y0)*100

    ##### Bias and relative bias for mu
    BIAS_mu[ii] <- mean(mu.0 - fit.y0)
    RBias_mu[ii] <- median(abs((mu.0 - fit.y0)/mu.0)*100)
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(ii)
}

result <- cbind(mean(BIAS_mu, na.rm = TRUE), 
                mean(RBias_mu, na.rm = TRUE), 
                mean(PI_y0, na.rm = TRUE))
result_df <- as.data.frame(result)
colnames(result_df) <- c("BIAS_mu",
                      "RBias_mu",
                      "PI_y0")
result_df

