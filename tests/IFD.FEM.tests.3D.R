##########################################
############## TEST SCRIPT ###############
############### for DEPTH ################
############### 3D domain ################
##########################################
# in order to use R under gdb debugger, write
# $ sudo R -d gdb

rm(list=ls())
graphics.off()
library(fdaPDE)

data(sphere3Ddata)
mesh<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)
plot(mesh)
nnodes = dim(mesh$nodes)[1]
FEMbasis = create.FEM.basis(mesh)

# create FEM functions
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- sin(2*pi*mesh$nodes[,1]) +  2 * sin(2*pi*mesh$nodes[,2]) +  sin(2*pi*mesh$nodes[,3])
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=0.05*abs(1 - (-1)))
  err <- 0 # rnorm(n=nnodes, mean=0, sd=abs(1 - (-1)))
  data[, i] <- func$coeff + c + err
  consts <- c(consts, c)
}

# con errori e smoothing
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- sin(2*pi*mesh$nodes[,1]) +  2 * sin(2*pi*mesh$nodes[,2]) +  sin(2*pi*mesh$nodes[,3])
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=0.05*abs(1 - (-1)))
  err <- rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
  solution <- smooth.FEM(observations = func$coeff + c + err, FEMbasis = FEMbasis, lambda = 0.5)
  data[, i] <- solution$fit.FEM$coeff
  consts <- c(consts, c)
}

data <- cbind(func$coeff, data)
data <- cbind(data, func$coeff + 5)

max<-max(data)
min<-min(data)

plot(FEM(data[,1], FEMbasis), max, min)#, FEMbasis$mesh$faces)

# define a weight function
weights <- function(eval_points)
{
  w <- eval_points[, 1]*eval_points[, 2]
  return(w)
}

w <- weights(FEMbasis$mesh$nodes)

# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", plot = TRUE)
sol$ifd
