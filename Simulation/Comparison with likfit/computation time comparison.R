rm(list = ls())
library(geoR)

#Source the functions
file_list <- list.files(path = "functions", full.names = TRUE)
for (file in file_list) {
  source(file)
}

# Covariance functions:  circular, exponential, Matern, and spherical
cov.method <- "circular"

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

w1 <- runif(n = n, min = 0, max = 1)
w2 <- runif(n = n, min = 0, max = 1)
gridw <- cbind(w1, w2)

set.seed(0519)
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


######################################################

library(microbenchmark)

benchmark_lowrank <- microbenchmark(
  model <- Krig.Gauss(data = data, 
                      cov.method = cov.method, 
                      knots = knots),
  times = 10L
)


benchmark_geoR <- microbenchmark(
  geoR.fit <- likfit(coords = data[,c("w1","w2")],
                     data = data$y,
                     ini.cov.pars=c(0.5, 0.5),
                     lik.met = "REML",
                     cov.model = cov.method),
  times = 10L
)

print(benchmark_lowrank, unit ="s")
print(benchmark_geoR, unit ="s")