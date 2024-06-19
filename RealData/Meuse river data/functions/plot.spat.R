plot.spat <- function(model, n1 = 100, n2 = 100, 
                      text_size = 1, legend_size = 1){
  
  packages <- c("spatstat.geom", "sp", "fields")
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      message(paste("Package", package, "is required but not installed."))
    }
  }
  text_size <- text_size
  legend_size <- legend_size
  coord <- model$coord
  K <- model$covs$K
  covs <- model$covs
  v_mode <- model$v_mode
  xi_estim <- model$xi_estim
  
  cvxhull <- spatstat.geom::convexhull.xy(coord)
  grid1 <- seq(min(coord[,1]), max(coord[,1]), length.out = n1)
  grid2 <- seq(min(coord[,2]), max(coord[,2]), length.out = n2)
  wgrid <- expand.grid(grid1, grid2)
  x.grid <- wgrid[,1]
  y.grid <- wgrid[,2]
  
  inpolygon <- sp::point.in.polygon(x.grid,y.grid,c(cvxhull$bdry[[1]]$x,cvxhull$bdry[[1]]$x[1]),c(cvxhull$bdry[[1]]$y,cvxhull$bdry[[1]]$y[1]))
  w1.0 <- x.grid[inpolygon==1]
  w2.0 <- y.grid[inpolygon==1]
  
  ##################
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
  Cspat <- cbind(w1.0, w2.0, Zw.0)
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(Cspat - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  
  xi_spat <- tail(xi_estim, K + 2)
  fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)
  
  ######################
  matrix_data <- matrix(inpolygon, nrow = n1, byrow = FALSE)
  one_positions <- which(matrix_data == 1, arr.ind = TRUE)
  matrix_data[one_positions] <- fit.w0
  matrix_data[matrix_data == 0] <- NA
  fields::image.plot(grid1, grid2, matrix_data,
             xlab = "w1",
             ylab = "w2",
             cex.axis = text_size,
             cex.lab = text_size,
             axis.args=list(cex.axis=legend_size))
}
