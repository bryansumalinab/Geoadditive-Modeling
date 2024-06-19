rm(list = ls())
library(geoR)

# Covariance functions:  circular, exponential, Matern, and spherical
cov.method <- "circular"

sill <- 0.5 # sill parameter
phi <- 0.15 # range parameter
sigma2 <- 0.10 # error variance
n <- 1000 # sample size

N <- 250 # number of simulations
BIAS_mu <- RBias_mu <- rep(NA, N)
PI_y0 <- rep(NA, N)
beta0 <- 3

set.seed(0519)
# Loop through each j
for (ii in 1:N) {
  tryCatch({
    w1 <- runif(n = n, min = 0, max = 1)
    w2 <- runif(n = n, min = 0, max = 1)
    gridw <- cbind(w1, w2)
    simdata <- geoR::grf(cov.pars = c(sill, phi),
                         grid = gridw, 
                         cov.model = cov.method)
    
    er <- rnorm(n, mean = 0, sd = sqrt(sigma2))
    ysp <- simdata$data
    y <- beta0 + ysp + er
    ymean <- beta0 + ysp
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
    
    geoR.fit <- likfit(coords = data.est[,c("w1","w2")],
                       data = data.est$y,
                       ini.cov.pars=c(0.5, 0.5),
                       lik.met = "REML",
                       cov.model = cov.method)
    coord0 <- cbind(data.pred$w1, data.pred$w2)
    pred_geoR <- krige.conv(coords = data.est[,c("w1","w2")], 
                            data = data.est$y, 
                            locations = coord0, 
                            krige = krige.control(obj.model = geoR.fit))
    
    y0 <- data.pred$y
    mu.0 <- data.pred$ymean
    ################## Prediction ################
    fit.y0 <- pred_geoR$predict
    y0.var.geoR <- pred_geoR$krige.var
    quantiles_geoR <- rbind(fit.y0 - 1.96*sqrt(y0.var.geoR),
                            fit.y0 + 1.96*sqrt(y0.var.geoR))
    pi_y0 <- rep(NA, length(y0))
    for (i in 1:length(y0)) {
      pi_y0[i] <- as.numeric(ifelse(y0[i] >= quantiles_geoR[1, i] & y0[i] <= quantiles_geoR[2, i], 1, 0))
    }
    
    ###### Prediction interval for y0
    PI_y0[ii] <- mean(pi_y0)*100

    ##### Bias and relative bias for mu
    BIAS_mu[ii] <- mean(mu.0 - fit.y0)
    RBias_mu[ii] <- median(abs((mu.0 - fit.y0)/mu.0)*100)
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
  cat("ITERATION", ii, "\n")
}

result <- cbind(mean(BIAS_mu, na.rm = TRUE), 
                mean(RBias_mu, na.rm = TRUE), 
                mean(PI_y0, na.rm = TRUE))
result_df <- as.data.frame(result)
colnames(result_df) <- c("BIAS_mu",
                         "RBias_mu",
                         "PI_y0")
result_df
