rm(list = ls())
file_list <- list.files("functions/", pattern = "\\.R$", full.names = TRUE)
source_files <- function(file_list) {
  for (file in file_list) {
    source(file)
  }
}
# Call the function to source all files
source_files(file_list)

data(meuse, package = "sp")
zinc <- log(meuse$zinc)
dist <- meuse$dist
coord <- meuse[,c("x","y")]/1000
K <- max(20, min(nrow(coord)/4, 150))
set.seed(0519)
knot.sp <- fields::cover.design(R = cbind(coord[,1], 
                                          coord[,2]), 
                                nd = K)
knots <- as.matrix(knot.sp$design)


###########################################################
####### (a) BIC comparison for different covariance ########
###########################################################
smooth_var <- c("dist", "elev")
X_smooth <- meuse[,smooth_var,drop = FALSE]
covar <- c("circular", 
           "exponential", 
           "matern", 
           "spherical")

BIC <- matrix(nrow = 4, ncol = 1)
rownames(BIC) <- covar
colnames(BIC) <- c("BIC")

set.seed(0519)
for (i in 1:length(covar)) {
  model <- Bsplinekrig.Gauss(y = zinc,
                             X_smooth = X_smooth,
                             X_linear = NULL,
                             coord = coord,
                             knots = knots,
                             cov.method = covar[i])
  BIC[i, 1] <- model$BIC
  cat("Index (i):", i, "\n")
  cat("Covariance type:", covar[i], "\n")
}

BIC


###########################################################
##### (b) Elevation as linear covariate using circular covariance
###########################################################
set.seed(0519)
model_elev_lin <- Bsplinekrig.Gauss(y = zinc,
                               X_smooth = meuse[,"dist",drop = FALSE],
                               X_linear = meuse[,"elev",drop = FALSE],
                               coord = coord,
                               knots = knots,
                               cov.method = "circular")
model_elev_lin$BIC


###########################################################
##### (c) Final model - using circular covariance
###########################################################
set.seed(0519)
model <- Bsplinekrig.Gauss(y = zinc,
                               X_smooth = meuse[,c("dist", "elev"),drop = FALSE],
                               X_linear = NULL,
                               coord = coord,
                               knots = knots,
                               cov.method = "circular")

# Estimated effect smooth term 1
plot.smooth(model, term = 1)
# Estimated effect smooth term 2
plot.smooth(model, term = 2)

# Fitted values
plot(model$fitted.values, zinc,
     xlab = "Fitted values",
     ylab = "Observed log(zinc)",
     main = "",
     pch = 19)
abline(a = 0, b = 1,
       col = "black", lwd = 2)

# Predicted spatial surface
plot.spat(model)

