# Supplementary material for:
# A roughness penalty approach to estimate densities over two-dimensional manifolds
# by Eleonora Arnone, Federico Ferraccioli, Clara Pigolotti, Laura M. Sangalli
# Computational Statistics and Data Analysis


# Plot a FEM object with jet colormap with range [m,M]
plot.FEM = function(FEM, M=NULL, m=NULL, ...){

  if (is.null(m)) { m = min(FEM$coeff[!is.na(FEM$coeff)])}
  if (is.null(M)) { M = max(FEM$coeff[!is.na(FEM$coeff)])}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order=FEM$FEMbasis$mesh$order
  nodes=FEM$FEMbasis$mesh$nodes
  edges=matrix(rep(0,6*ntriangles),ncol=2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,]=c(triangles[3*order*i+1],triangles[3*order*i+2])
    edges[3*i+2,]=c(triangles[3*order*i+1],triangles[3*order*i+3])
    edges[3*i+3,]=c(triangles[3*order*i+2],triangles[3*order*i+3])
  }
  edges=edges[!duplicated(edges),]
  edges<-as.vector(t(edges))

  coeff = FEM$coeff

  FEMbasis = FEM$FEMbasis

  mesh = FEMbasis$mesh

  p=jet.col(n=128,alpha=0.8)
  # p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(128)
  palette(p)

  ncolor=length(p)

  nplots <- ncol(coeff)

  for (i in 1:nplots) {

    if (i > 1)
      readline("Press any key for the next plot...")

    diffrange = M - m

    col = coeff[triangles,i]
    col = (col - m)/diffrange*(ncolor-1)+1
    col[is.na(col)] = "grey"

    open3d()
    axes3d()

    z <- coeff[triangles,i]
    rgl.triangles(nodes[triangles,1], nodes[triangles,2], z,
                  color = col,...)

    aspect3d(2,2,1)
    rgl.viewpoint(0,-45)
  }
}

zoom = 1
userMatrix = rbind(c( 0.96563137,  0.1774523, -0.1899119,    0),
                   c( -0.03294301,  0.8083354,  0.5877997,    0),
                   c( 0.25781897, -0.5613416,  0.7863999,    0),
                   c(0.00000000,  0.0000000,  0.0000000 ,   1))
windowRect = c(150,  150, 420,  420)
