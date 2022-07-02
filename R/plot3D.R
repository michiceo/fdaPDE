# Supplementary material for:
# A roughness penalty approach to estimate densities over three-dimensional manifolds
# by Eleonora Arnone, Federico Ferraccioli, Clara Pigolotti, Laura M. Sangalli
# Computational Statistics and Data Analysis


# Plot a FEM object with jet colormap with range [m,M]
plot.FEM = function(FEM, M=NULL, m=NULL, ...){

  if (is.null(m)) { m = min(FEM$coeff[!is.na(FEM$coeff)])}
  if (is.null(M)) { M = max(FEM$coeff[!is.na(FEM$coeff)])}
  order=FEM$FEMbasis$mesh$order
  nodes=FEM$FEMbasis$mesh$nodes
  if (class(FEM$FEMbasis$mesh) == "mesh.3D"){
    triangles = c(t(FEM$FEMbasis$mesh$tetrahedrons))
    ntriangles = nrow(FEM$FEMbasis$mesh$tetrahedrons)

    faces=matrix(rep(0,12*ntriangles),ncol=3)
    for(i in 0:(ntriangles-1)){
      faces[4*i+1,]=c(triangles[4*order*i+1],triangles[4*order*i+2],triangles[4*order*i+3])
      faces[4*i+2,]=c(triangles[4*order*i+1],triangles[4*order*i+2],triangles[4*order*i+4])
      faces[4*i+3,]=c(triangles[4*order*i+1],triangles[4*order*i+3],triangles[4*order*i+4])
      faces[4*i+4,]=c(triangles[4*order*i+2],triangles[4*order*i+3],triangles[4*order*i+4])
    }
    faces=faces[!duplicated(faces),]
    edges<-as.vector(t(faces))

    #edges <- as.vector(t(FEMbasis$mesh$faces[as.logical(FEMbasis$mesh$facesmarkers),]))
  }
  else{
    triangles = c(t(FEM$FEMbasis$mesh$triangles))
    ntriangles = nrow(FEM$FEMbasis$mesh$triangles)

    edges=matrix(rep(0,6*ntriangles),ncol=2)
    for(i in 0:(ntriangles-1)){
      edges[3*i+1,]=c(triangles[3*order*i+1],triangles[3*order*i+2])
      edges[3*i+2,]=c(triangles[3*order*i+1],triangles[3*order*i+3])
      edges[3*i+3,]=c(triangles[3*order*i+2],triangles[3*order*i+3])
    }
    edges=edges[!duplicated(edges),]
    edges<-as.vector(t(edges))
  }

  coeff = FEM$coeff

  FEMbasis = FEM$FEMbasis

  mesh = FEMbasis$mesh

  p=jet.col(n=128,alpha=0.8)
  # p <- colorRampPalette(c("#0E1E44","#3E6DD8","#68D061","#ECAF53", "#EB5F5F","#E11F1C"))(128)
  palette(p)

  ncolor=length(p)

  if (class(FEM$FEMbasis$mesh) == "mesh.3D"){
    nsurf = dim(coeff)[[2]]
    for (isurf in 1:nsurf)
    {
      open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
      rgl.pop("lights")
      light3d(specular="black")

      diffrange = M - m

      col = coeff[edges,isurf]
      col = (col - min(coeff[,isurf][!is.na(coeff[,isurf])]))/diffrange*(ncolor-1)+1
      col[is.na(col)] = "grey"

      rgl.triangles(x = nodes[edges ,1], y = nodes[edges ,2],
                    z=nodes[edges,3],
                    color = col,...)
      # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
      #           z=nodes[edges,3],
      #           color = "black",...)
      aspect3d("iso")

      if (nsurf > 1 && isurf<nsurf)
      {readline("Press a button for the next plot...")}
    }
  }
  else{
    nsurf = dim(coeff)[[2]]
    for (isurf in 1:nsurf)
    {
      open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
      rgl.pop("lights")
      light3d(specular="black")

      diffrange = M - m

      col = coeff[triangles,isurf]
      col = (col - min(coeff[,isurf][!is.na(coeff[,isurf])]))/diffrange*(ncolor-1)+1
      col[is.na(col)] = "grey"

      rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                    z=nodes[triangles,3],
                    color = col,...)
      # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
      #           z=nodes[edges,3],
      #           color = "black",...)
      aspect3d("iso")

      if (nsurf > 1 && isurf<nsurf)
      {readline("Press a button for the next plot...")}
    }
  }
}

zoom = 1
userMatrix = rbind(c( 0.96563137,  0.1774523, -0.1899119,    0),
                   c( -0.03294301,  0.8083354,  0.5877997,    0),
                   c( 0.25781897, -0.5613416,  0.7863999,    0),
                   c(0.00000000,  0.0000000,  0.0000000 ,   1))
windowRect = c(150,  150, 420,  420)
