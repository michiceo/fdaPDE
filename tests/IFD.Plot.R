#if(2D)
  
plot.image.2D = function(fun, max, min){
  if(class(fun$FEMbasis$mesh)!='mesh.2D')
    stop('This function is implemented only for 2D mesh FEM objects')
  n = length(fun$FEMbasis$mesh$nodes[!duplicated(fun$FEMbasis$mesh$nodes[,1]),1]);
  m = length(fun$FEMbasis$mesh$nodes[!duplicated(fun$FEMbasis$mesh$nodes[,2]),2]);
  x = sort(fun$FEMbasis$mesh$nodes[!duplicated(fun$FEMbasis$mesh$nodes[,1]),1]);
  y = sort(fun$FEMbasis$mesh$nodes[!duplicated(fun$FEMbasis$mesh$nodes[,2]),2]);
  location = expand.grid(x,y);
  z = matrix(eval.FEM(fun,location),nrow=n,ncol=m);
 
  # xx = (x - min(x))/(max(x) - min(x));
  # yy = (y - min(y))/(max(y) - min(y));
  zlim = c(min, max)

  image(x,y,z,zlim, axes = FALSE)
  #contour(x,y,z, add = TRUE)
}

# plot image with contour of the function VS lower whisker / Q1
plot.image.diff_lw_q1.2D = function(func,lw_q1, max, min){
  if(class(func$FEMbasis$mesh)!='mesh.2D')
    stop('This function is implemented only for 2D mesh FEM objects')
  if(class(lw_q1$FEMbasis$mesh)!='mesh.2D')
    stop('This function is implemented only for 2D mesh FEM objects')
  
  n = length(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,1]),1]);
  m = length(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,2]),2]);
  x = sort(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,1]),1]);
  y = sort(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,2]),2]);
  location = expand.grid(x,y);
  zz = matrix(eval.FEM(func,location) - eval.FEM(lw_q1,location),nrow=n,ncol=m);
  z = matrix(eval.FEM(func,location),nrow=n,ncol=m);
  z[which(zz < 0)] = NA;
  zlim = c(min, max)
  
  image(x,y,z,zlim, axes = FALSE)
  #contour(x,y,z, add = TRUE)
}

# plot image with contour of the function VS upper whisker / Q3
plot.image.diff_uw_q3.2D = function(func,uw_q3, max, min){
  if(class(func$FEMbasis$mesh)!='mesh.2D')
    stop('This function is implemented only for 2D mesh FEM objects')
  if(class(uw_q3$FEMbasis$mesh)!='mesh.2D')
    stop('This function is implemented only for 2D mesh FEM objects')
  
  n = length(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,1]),1]);
  m = length(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,2]),2]);
  x = sort(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,1]),1]);
  y = sort(func$FEMbasis$mesh$nodes[!duplicated(func$FEMbasis$mesh$nodes[,2]),2]);
  location = expand.grid(x,y);
  zz = matrix(eval.FEM(uw_q3,location) - eval.FEM(func,location),nrow=n,ncol=m);
  z = matrix(eval.FEM(func,location),nrow=n,ncol=m);
  z[which(zz < 0)] = NA;
  zlim = c(min, max)
  
  image(x,y,z,zlim, axes = FALSE)
  #contour(x,y,z, add = TRUE)
}

#if(3D)

plot.diff_lw_q1.3D = function(func, lw_q1, max, min){
  if(class(func$FEMbasis$mesh)!='mesh.2.5D' & class(func$FEMbasis$mesh)!='mesh.3D')
    stop('This function is implemented only for 2.5D/3D mesh FEM objects')
  
  diff = func$coeff - lw_q1$coeff
  func$coeff[which(diff < 0)] = NA
  
  plot(func, max, min)
}

plot.diff_uw_q3.3D = function(func, uw_q3, max, min){
  if(class(func$FEMbasis$mesh)!='mesh.2.5D' & class(func$FEMbasis$mesh)!='mesh.3D')
    stop('This function is implemented only for 2.5D/3D mesh FEM objects')
  diff = uw_q3$coeff - func$coeff
  func$coeff[which(diff < 0)] = NA
  
  plot(func, max, min)
}