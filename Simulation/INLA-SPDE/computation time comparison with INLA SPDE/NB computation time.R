rm(list = ls())

#Source the functions
file_list <- list.files(path = "functionsNB", full.names = TRUE)
for (file in file_list) {
  source(file)
}

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

######### Computation time ##############
library(microbenchmark)

##### Computation time for LPS######
benchmark_LPS <- microbenchmark(
  {################################################
    knot.sp <- fields::cover.design(R = cbind(runif(n = n, min = minysp, max = maxysp), 
                                              runif(n = n, min = minysp, max = maxysp)), 
                                    nd = 150)
    knots <- as.matrix(knot.sp$design)
    
    
    model <- Krig.NB(data = dat,
                     cov.method = "exponential",
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
    
    ########## New data for prediction ###################
    
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
    
    ###### smooth component
    ngrid <- 100
    xgrid <- seq(min(xsmooth), max(xsmooth), length.out = ngrid)
    ind <- (3): (2 + nbasis)
    Zsm.grid <- fsmooth.pred.bspline(x0 = xgrid,
                                     B = nbasis,
                                     const = const.sm,
                                     xmin = min(xsmooth),
                                     xmax = max(xsmooth))
    
    xlin.0 <- runif(ngrid)
    C0 <- cbind(1, xlin.0, Zsm.grid, Cspat.cent)
    
    ########### Posterior predictive samples
    # Compute mean and variance for log(mu)
    logmu <- as.numeric(C0%*%xi_estim)
    logmu.var <- as.numeric(diag(C0%*%Sigma%*%t(C0)))
    
    # Generate NB samples
    r.NB <- list()
    N <- 1000
    for (i in 1:length(logmu)) {
      rn <- rnorm(N, mean = logmu[i], sd = sqrt(logmu.var[i]))
      mu <- exp(rn)
      r.NB[[i]] <- MASS::rnegbin(n = N, mu = mu, theta = exp(v_mode[length(v_mode) - 1]))
    }},
  times = 10L
)

################ INLA-SPDE ###########################
library(INLA)
library(fmesher)
coord <- data.frame(w1, w2)
y <- dat$y

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

##### Computation time for INLA-SPDE######
benchmark_SPDE <- microbenchmark(
  {### spatial component
    
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
    
    #### predictive samples
    index_pred.ymean <- inla.stack.index(stack, "pred.ymean")$data
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
    }},
  times = 10L
)

# results
summary(benchmark_LPS, unit = "s")[, -1]
summary(benchmark_SPDE, unit = "s")[, -1]
