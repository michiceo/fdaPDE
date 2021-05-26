CPP_FEM.IFD <- function(data, FEMbasis, ndim, mydim, weights, search)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(weights) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(search) <- "integer"

  ## Call C++ function
  bigsol <- .Call("Integrated_Functional_Depth", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, weights, search,
                  PACKAGE = "fdaPDE")

  ## Reset them correctly
  # FEMbasis$mesh$triangles = FEMbasis$mesh$triangles + 1
  # FEMbasis$mesh$edges = FEMbasis$mesh$edges + 1
  # FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] + 1
  
  return(bigsol)
}


CPP_FEM.manifold.IFD <- function(data, FEMbasis, ndim, mydim, weights, search)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  
  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(weights) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(search) <- "integer"
  
  ## Call C++ function
  bigsol <- .Call("Integrated_Functional_Depth", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, weights, search,
                  PACKAGE = "fdaPDE")
  
  ## Reset them correctly
  # FEMbasis$mesh$triangles = FEMbasis$mesh$triangles + 1
  # FEMbasis$mesh$edges = FEMbasis$mesh$edges + 1
  # FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] + 1
  
  return(bigsol)
}


CPP_FEM.volume.IFD <- function(data, FEMbasis, ndim, mydim, weights, search)
{
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  ## Set proper type for correct C++ reading
  storage.mode(data) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(weights) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(search) <- "integer"
  
  bigsol <- .Call("Integrated_Functional_Depth", data, FEMbasis$mesh, FEMbasis$order, mydim, ndim, weights, search,
                  PACKAGE = "fdaPDE")
  
  ## Reset them correctly
  # FEMbasis$mesh$triangles = FEMbasis$mesh$triangles + 1
  # FEMbasis$mesh$edges = FEMbasis$mesh$edges + 1
  # FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] + 1
  
  return(bigsol)
}

