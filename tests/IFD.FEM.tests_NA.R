##########################################
############## TEST SCRIPT ###############
############### for DEPTH ################
##########################################

rm(list=ls())
graphics.off()
library(fdaPDE)
source("censoring.R")

set.seed(5847947)

a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)

############################## 2D ##############################

#### Test 1: 2D square domain ####
#            locations = nodes
#            order FE = 1
library(fdaPDE)
source("plot2D.R")

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

### 1. create FEM functions
func1 = FEM(0.9*mesh$nodes[, 1], FEMbasis)
func2 = FEM(1.1*mesh$nodes[, 1], FEMbasis)
func3 = FEM(1.5*mesh$nodes[, 1], FEMbasis)
# add an error component to the functions
set.seed(5847947)
truedatarange<-max(c(func1$coeff,func2$coeff,func3$coeff))-min(c(func1$coeff,func2$coeff,func3$coeff))
err1<-rnorm(n=nnodes, mean=0, sd=0.005)
err2<-rnorm(n=nnodes, mean=0, sd=0.005)
err3<-rnorm(n=nnodes, mean=0, sd=0.005)
func1$coeff <- func1$coeff + err1
func2$coeff <- func2$coeff + err2
func3$coeff <- func3$coeff + err3

# build data matrix
data <- cbind(func1$coeff, func2$coeff, func3$coeff)
# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", nThreads = 5)
ifd <- sol$ifd

max<-max(data[!is.na(data)])
min<-min(data[!is.na(data)])

# create missing observations in the data
# data <- createNA(Data = data, schema = 2, mesh_ref = mesh, locations = locations)
# schema = 2 tanto laborioso per q=100
ifd_NA <- ifd
num_NA <- 0
for(it in 1:5){
  data <- createNA(Data = data, p=0.99, schema = 4, mesh_ref = mesh, locations = locations) # p=0.75 meno NA -> p=0.5 più NA
  
  max<-max(data[!is.na(data)])
  min<-min(data[!is.na(data)])
  
  # compute the solution to the problem
  # (compute the depth value for the set of functions considered)
  sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", nThreads = 5)
  ifd_NA <- cbind(ifd_NA, sol$ifd)
  num_NA <- cbind(num_NA, sum(is.na(data)))
}

plot(num_NA[1,], ifd_NA[1,], main = "NA number VS ranking", xlab = "NA number(num_NA)", ylab = "depth value for function i (ifd[i])", col = "red", type = "l")
lines(num_NA[1,], ifd_NA[2,], col = "green")
lines(num_NA[1,], ifd_NA[3,], col = "blue")
legend(x = "topright", legend = c("f(x) = x", "f(x) = 0.9-x", "f(x) = 0.1"), lty = c(1,1,1), col = c("black","red","green") )

plot.image.2D(FEM(data[,1], FEMbasis), max, min)
plot.image.2D(FEM(data[,2], FEMbasis), max, min)
plot.image.2D(FEM(data[,3], FEMbasis), max, min)



### 2. generate a bundle of functions
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
data <- cbind(data, func$coeff + 5)

# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", plot=TRUE)
ifd <- sol$ifd

max<-max(data[!is.na(data)])
min<-min(data[!is.na(data)])

plot.image.2D(FEM(data[,1], FEMbasis), max, min)

# create missing observations in the data
# data <- createNA(Data = data, schema = 2, mesh_ref = mesh, locations = locations)
# schema = 2 tanto laborioso per q=100
ifd_NA <- ifd
num_NA <- 0
for(it in 1:5){
  data <- createNA(Data = data, p=0.99, schema = 4, mesh_ref = mesh, locations = locations) # p=0.75 meno NA -> p=0.5 più NA
  
  max<-max(data[!is.na(data)])
  min<-min(data[!is.na(data)])
  
  # compute the solution to the problem
  # (compute the depth value for the set of functions considered)
  sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD", plot=TRUE)
  ifd_NA <- cbind(ifd_NA, sol$ifd)
  num_NA <- cbind(num_NA, sum(is.na(data)))
}

plot(num_NA[1,], ifd_NA[1,], main = "NA number VS ranking", xlab = "NA number(num_NA)", ylab = "depth value for function i (ifd[i])", col = "red", type = "l")
lines(num_NA[1,], ifd_NA[2,], col = "green")
lines(num_NA[1,], ifd_NA[3,], col = "blue")
legend(x = "topright", legend = c("f(x) = x", "f(x) = 0.9-x", "f(x) = 0.1"), lty = c(1,1,1), col = c("black","red","green") )

plot.image.2D(FEM(data[,1], FEMbasis), max, min)

############################## 2.5D ##############################

#### Test 3: 2.5D horseshoe domain ####
#            locations = nodes
#            order FE = 1
#rm(list=ls())
#graphics.off()
library(fdaPDE)

# generate the mesh
data(horseshoe2.5D)
mesh<-create.mesh.2.5D(horseshoe2.5D$nodes, horseshoe2.5D$triangles)
plot(mesh)
nnodes = dim(mesh$nodes)[1]
FEMbasis = create.FEM.basis(mesh)

# create FEM functions
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- fs.test.3D(mesh$nodes[,1],mesh$nodes[,2],mesh$nodes[,3])
func <- FEM(f, FEMbasis)
for(i in 1:nn){
  c <- rnorm(1, mean=0, sd=2)
  err <- 0 # rnorm(n=nnodes, mean=0, sd=abs(1 - (-1)))
  data[, i] <- func$coeff + c + err
  consts <- c(consts, c)
}

# con errori e smoothing
nn <- 50
data <- matrix(0, nnodes, nn)
consts <- NULL
f <- fs.test.3D(mesh$nodes[,1],mesh$nodes[,2],mesh$nodes[,3])
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

# compute the solution to the problem
# (compute the depth value for the set of functions considered)
sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD")
ifd <- sol$ifd # sembrano molto piccoli i numeri, ma vengono piccoli sia senza che con NA

max<-max(data[!is.na(data)])
min<-min(data[!is.na(data)])

plot(FEM(data[,1], FEMbasis), max, min)

# create missing observations in the data
# data <- createNA(Data = data, schema = 2, mesh_ref = mesh, locations = locations)
# schema = 2 tanto laborioso per q=100
ifd_NA <- ifd
num_NA <- 0
for(it in 1:5){
  data <- createNA(Data = data, p=0.99, schema = 4, mesh_ref = mesh, locations = locations) # p=0.75 meno NA -> p=0.5 più NA

  max<-max(data[!is.na(data)])
  min<-min(data[!is.na(data)])

  plot(FEM(data[,1], FEMbasis), max, min)

  # compute the solution to the problem
  # (compute the depth value for the set of functions considered)
  sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MBD")
  ifd_NA <- cbind(ifd_NA, sol$ifd)
  num_NA <- cbind(num_NA, sum(is.na(data)))
}

plot(num_NA[1,], ifd_NA[1,], main = "NA number VS ranking", xlab = "NA number(num_NA)", ylab = "depth value for function i (ifd[i])", col = "red", type = "l")
lines(num_NA[1,], ifd_NA[2,], col = "green")
lines(num_NA[1,], ifd_NA[3,], col = "blue")
legend(x = "topright", legend = c("f(x) = x", "f(x) = 0.9-x", "f(x) = 0.1"), lty = c(1,1,1), col = c("black","red","green") )

plot(FEM(data[,1], FEMbasis), max, min)
