##########################################
############## TEST SCRIPT ###############
############### for DEPTH ################
############### 2D domain ################
##########################################
# in order to use R under gdb debugger, write
# $ sudo R -d gdb

rm(list=ls())
graphics.off()
library(fdaPDE)

set.seed(5847947)

a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)

#### Test 1: 2D square domain ####
#            locations = nodes
#            order FE = 1

# generate the mesh
q <- 100
x = seq(0,1, length.out = q)
y = x
locations = expand.grid(x,y)
mesh = create.mesh.2D(locations) #, order=2)
plot(mesh)
nnodes = dim(mesh$nodes)[1]
FEMbasis = create.FEM.basis(mesh)

# create censoring with censoring.R for the following functions

##### 1. nn = 3 linear functions + err #####
# create FEM functions
func1 = FEM(0.9*mesh$nodes[, 1], FEMbasis) # FEM(mesh$nodes[, 1], FEMbasis)
func2 = FEM(1.1*mesh$nodes[, 1], FEMbasis) # FEM(0.9-mesh$nodes[, 1], FEMbasis)
func3 = FEM(1.5*mesh$nodes[, 1], FEMbasis) # FEM(rep(0.1, nnodes), FEMbasis)
# add an error component to the functions
set.seed(5847947)
truedatarange<-max(c(func1$coeff,func2$coeff,func3$coeff))-min(c(func1$coeff,func2$coeff,func3$coeff))
err1<-rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
err2<-rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
err3<-rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
func1$coeff <- func1$coeff + err1
func2$coeff <- func2$coeff + err2
func3$coeff <- func3$coeff + err3

# build data matrix
data <- cbind(func1$coeff, func2$coeff, func3$coeff)

max<-max(data)
min<-min(data)

plot.image.2D(func1, max, min)
plot.image.2D(func2, max, min)
plot.image.2D(func3, max, min)

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



##### 2. bundle of sinusoidal functions + err #####
# create FEM functions
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- a1*sin(2*pi*mesh$nodes[,1])*cos(2*pi*mesh$nodes[,2])+a2*sin(3*pi*mesh$nodes[,1])
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=2)
  err <- 0 # rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
  data[, i] <- func$coeff + c + err
  consts <- c(consts, c)
}

# con errori e smoothing
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- a1*sin(2*pi*mesh$nodes[,1])*cos(2*pi*mesh$nodes[,2])+a2*sin(3*pi*mesh$nodes[,1])
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=0.05*abs(1 - (-1)))
  err <- rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
  solution <- smooth.FEM(observations = func$coeff + c + err, FEMbasis = FEMbasis, lambda = 0.5)
  data[, i] <- solution$fit.FEM$coeff
  consts <- c(consts, c)
}

data <- cbind(func$coeff, data)

data <- cbind(data, func$coeff+5)

max<-max(data)
min<-min(data)

plot.image.2D(FEM(data[,1], FEMbasis), max, min)

# use perspective plot
plotcolors <- palette(rainbow(51))
persp3d(x, y, f, theta = 30, phi = 30, expand = 0.5, shade = 0.1,
        col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f")
for(i in 2:51){
  persp3d(x, y, func$coeff + consts[i], theta = 30, phi = 30, expand = 0.5, shade = 0.1,
          col = plotcolors[i], xlab = "x", ylab = "y", zlab = "f", add=T)
}

# plot just 5 of them
plotcolors <- c("gray", palette(rainbow(51)))
persp3d(x, y, f, theta = 30, phi = 30, expand = 0.5, shade = 0.1,
        col = plotcolors[2], xlab = "x", ylab = "y", zlab = "f")
for(i in 2:5){
  persp3d(x, y, data[, i], theta = 30, phi = 30, expand = 0.5, shade = 0,
          col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f", add=T)
}

# define a weight function
weights <- function(eval_points)
{
  w <- eval_points[, 1]*eval_points[, 2]
  return(w)
}

w <- weights(FEMbasis$mesh$nodes)

# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", plot=TRUE, func = data[,51])
sol$ifd

# plot the functional boxplot in the "2D" way (not possible in higher dimensions)
plotcolors <- palette(rainbow(5))
persp3d(x, y, sol$median, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f(x)")
persp3d(x, y, sol$firstQuartile, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[2], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$thirdQuartile, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[3], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$lowerWhisker, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[4], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$upperWhisker, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[5], xlab = "x", ylab = "y", zlab = "f(x)", add=T)

legend3d("topright", legend = c('median', 'first quartile', 'third quartile', 'lower whisker', 'upper whisker'), pch = 16, col = rainbow(5), cex=3, inset=c(0.1))
setwd("/home/michiceo/Immagini/fdaPDE_simulations")
snapshot3d(filename = 'fb_2D.png', fmt = 'pdf')



##### 3. bundle of gaussian functions + err #####
# create FEM functions
library(mvtnorm)
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
s <- diag(rep(0.01, 2), 2, 2, names=TRUE)
f <- dmvnorm(mesh$nodes, mean = rep(0, 2), sigma = s)
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=4)
  err <- rnorm(n=nnodes, mean=0, sd=0.05*abs(1 - (-1)))
  data[, i] <- func$coeff + c + err
  consts <- c(consts, c)
}
data <- cbind(func$coeff, data)

max<-max(data)
min<-min(data)

plot.image.2D(FEM(data[,1], FEMbasis), max, min)

# use perspective plot
plotcolors <- palette(rainbow(51))
persp3d(x, y, f, theta = 30, phi = 30, expand = 0.5, shade = 0.1,
        col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f")
for(i in 2:51){
  persp3d(x, y, func$coeff + consts[i], theta = 30, phi = 30, expand = 0.5, shade = 0.1,
          col = plotcolors[i], xlab = "x", ylab = "y", zlab = "f", add=T)
}

plotcolors <- palette(rainbow(51))
persp3d(x, y, f, theta = 30, phi = 30, expand = 0.5, shade = 0.1,
        col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f")
for(i in 2:51){
  persp3d(x, y, data[, i], theta = 30, phi = 30, expand = 0.5, shade = 0.1,
          col = plotcolors[i], xlab = "x", ylab = "y", zlab = "f", add=T)
}

# define a weight function
weights <- function(eval_points)
{
  w <- eval_points[, 1]*eval_points[, 2]
  return(w)
}

w <- weights(FEMbasis$mesh$nodes)

# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, weights, depth_choice = "MBD")
sol$ifd

# plot the functional boxplot in the "2D" way (not possible in higher dimensions)
plotcolors <- palette(rainbow(5))
persp3d(x, y, sol$median, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[1], xlab = "x", ylab = "y", zlab = "f(x)")
persp3d(x, y, sol$firstQuartile, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[2], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$thirdQuartile, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[3], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$lowerWhisker, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[4], xlab = "x", ylab = "y", zlab = "f(x)", add=T)
persp3d(x, y, sol$upperWhisker, theta = 60, phi = 30, expand = 3, shade = 0.1,
        col = plotcolors[5], xlab = "x", ylab = "y", zlab = "f(x)", add=T)

legend3d("topright", legend = c('median', 'first quartile', 'third quartile', 'lower whisker', 'upper whisker'), pch = 16, col = rainbow(5), cex=3, inset=c(0.1))
setwd("/home/michiceo/Immagini/fdaPDE_simulations")
snapshot3d(filename = 'fb_2D.png', fmt = 'pdf')
