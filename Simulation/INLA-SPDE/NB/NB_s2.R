rm(list = ls())
library(INLA)
library(fmesher)

# Function s2
smooth_spat <- function(x, y) {
  return((x^3 + x*y + y^2)/25)
}

minysp <- -3
maxysp <- 3
n <- 1000
sderr <- 0.25
beta0 <- 3
beta1 <- -0.5

N <- 250

CI_sm <- RBias_sm <- numeric(N)
PI_y0 <- numeric(N)
Bias_mu <- RBias_mu <- numeric(N)

set.seed(0519)

for (ii in 1:N) {
  tryCatch({
    ####### Data for estimation
    xsmooth <- runif(n, -1, 1)
    fsm <- function(x) cos(2*pi*x)
    fsmooth <- fsm(xsmooth)
    
    w1 <- runif(n, min = minysp, max = maxysp)
    w2 <- runif(n, min = minysp, max = maxysp)
    xlin <- runif(n)
    ysp <- smooth_spat(x = w1, y = w2)
    er <- rnorm(n = n, mean = 0, sd = sderr)
    ymean <- exp(beta0 + beta1 * xlin + fsmooth + ysp)
    y <- rpois(n = n, 
               lambda = ymean * exp(er))
    coord <- data.frame(w1, w2)
    
    
    ####### Data for prediction
    x1 <- seq(minysp, maxysp, length.out = 10)
    x2 <- seq(minysp , maxysp, length.out = 10)
    grid <- expand.grid(x1, x2)
    w1.0 <- grid[,1]
    w2.0 <- grid[,2]
    coord.0 <- data.frame("w1" = w1.0, "w2" = w2.0)
    ysp.0 <- smooth_spat(x = w1.0, y = w2.0)
    ngrid <- length(ysp.0)
    xgrid <- seq(min(xsmooth), max(xsmooth), length.out = ngrid)
    xlin.0 <- runif(length(ysp.0))
    er.0 <- rnorm(n = ngrid, mean = 0, sd = sderr)
    ymean.0 <- exp(beta0 + beta1 * xlin.0 + fsm(xgrid) + ysp.0)
    y0 <- rpois(n = n, 
                lambda = ymean.0 * exp(er.0))
    
    ### spatial component
    
    # (1) create mesh
    mesh <- fm_mesh_2d_inla(loc = coord,
                            max.edge = 1) 
    spde <- inla.spde2.matern(mesh)
    
    # (2) Basis matrix
    A.est <- inla.spde.make.A(mesh = mesh,
                              loc = as.matrix(coord))
    # (3) For vector of coefficients
    f.index <- inla.spde.make.index(name = "s.field",
                                    n.spde = spde$n.spde)
    
    #### smooth component
    # (1) create mesh
    mesh1 <- fm_mesh_1d(xsmooth)
    spde1 <- inla.spde2.matern(mesh1, constr = TRUE)
    
    # (2) Basis matrix
    A.est1 <- inla.spde.make.A(mesh = mesh1,
                               loc = xsmooth)
    # (3) For vector of coefficients
    f.index1 <- inla.spde.make.index(name = "xs",
                                     n.spde = spde1$n.spde)
    
    stack.est <- inla.stack(
      data = list(y = y),
      A = list(A.est, A.est1, 1, 1),
      effects = list(f.index, f.index1, 
                     list(Intercept = rep(1, dim(A.est)[1])),
                     xlin = xlin),
      tag = "est")
    
    A.pred <- inla.spde.make.A(mesh = mesh,
                               loc = as.matrix(coord.0))
    A.pred1 <- inla.spde.make.A(mesh = mesh1,
                                loc = xgrid)
    
    # mean model pred
    stack.pred.ymean <- inla.stack(
      data = list(y = NA),
      A = list(A.pred, A.pred1, 1, 1),
      effects = list(f.index, f.index1, 
                     list(Intercept = rep(1, dim(A.pred)[1])),
                     list(xlin = xlin.0)),
      tag = "pred.ymean")
    
    # smooth component
    stack.pred.smooth <- inla.stack(
      data = list(y = NA),
      A = list(A.pred1),
      effects = list(f.index1),
      tag = "pred.smooth")
    
    
    # Combine the stacks for fitting and predictions
    stack <- inla.stack(stack.est, 
                        stack.pred.ymean, 
                        stack.pred.smooth)
    
    formula <- y ~ -1 + Intercept + xlin + f(xs, model = spde1) + f(s.field, model = spde)
    
    
    inla.mod <- inla(formula, 
                     data = inla.stack.data(stack, spde = spde),
                     family="nbinomial",
                     control.predictor = list(A = inla.stack.A(stack), 
                                              compute = TRUE, link = 1 ),
                     control.compute=list(config = TRUE)) # required for inla.posterior.sample()
    
    
    ########## Prediction
    
    ##### ymean prediction
    index_pred.ymean <- inla.stack.index(stack, "pred.ymean")$data
    fitted.ymean   <- inla.mod$summary.fitted.values$mean[index_pred.ymean]
    fitted.ymean.lowerci <- inla.mod$summary.fitted.values$`0.025quant`[index_pred.ymean]
    fitted.ymean.upperci <- inla.mod$summary.fitted.values$`0.975quant`[index_pred.ymean]
    
    
    ##### smooth prediction
    index_pred.smooth <- inla.stack.index(stack, "pred.smooth")$data
    fitted.smooth   <- inla.mod$summary.linear.predictor$mean[index_pred.smooth]
    fitted.smooth.lowerci <- inla.mod$summary.linear.predictor$`0.025quant`[index_pred.smooth]
    fitted.smooth.upperci <- inla.mod$summary.linear.predictor$`0.975quant`[index_pred.smooth]
    
    
    #### smooth component
    fsxgrid <- fsm(xgrid) - mean(fsm(xgrid))
    ci_sm <- as.numeric(fsxgrid >= fitted.smooth.lowerci & fsxgrid <= fitted.smooth.upperci)
    CI_sm[ii] <- mean(ci_sm)*100
    RBias_sm[ii] <- median(abs((fsxgrid - fitted.smooth)/fsxgrid)*100)
    
    ### mean
    Bias_mu[ii] <- mean(ymean.0 - fitted.ymean)
    RBias_mu[ii] <- mean(abs((ymean.0 - fitted.ymean)/ymean.0)*100)
    
    ## predictive samples
    n.samples <- 500
    posterior.samples <- inla.posterior.sample(n.samples, inla.mod)
    pred.samples <- matrix(NA, nrow = n.samples, ncol = length(index_pred.ymean))
    for (i in 1:n.samples) {
      eta <- posterior.samples[[i]]$latent[index_pred.ymean]
      lambdapost <- exp(eta)
      thetapost <- as.numeric(posterior.samples[[i]]$hyperpar[1])
      pred.samples[i, ] <- rnbinom(length(lambdapost), 
                                   size = thetapost, 
                                   mu = lambdapost)
    }
    
    lower.y0 <- apply(pred.samples, 2, quantile, probs = 0.025)  # Lower 2.5% quantile
    upper.y0 <- apply(pred.samples, 2, quantile, probs = 0.975)  # Upper 97.5% quantile
    PI_y0[ii] <- mean(y0 >= lower.y0 & y0 <= upper.y0) * 100
    
  },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  print(ii)
}

# Extreme values are included in the results
result <- data.frame(
  Metric = c("Bias", "Relative Bias", "Coverage"),
  smooth = c(NA, mean(RBias_sm), mean(CI_sm)),
  mu = c(mean(Bias_mu[which(Bias_mu > -100)]), mean(RBias_mu[which(RBias_mu < 100)]), NA),
  y0 = c(NA, NA, mean(PI_y0))
)

result
