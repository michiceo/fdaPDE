#' Integrated functional depth computation for complex multidimensional functional data
#' @param data A matrix of dimensions #mesh nodes-by-#functions. Data are functions:
#' each row corresponds to the evaluations of the functions at one specific node of the mesh, 
#' each column corresponds to the evaluation of a specific function at all the mesh nodes.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis,
#' as created by \code{\link{create.FEM.basis}}.
#' @param weights A vector of length #\code{nodes} of the mesh. It corresponds to the
#' weights for the integration. The sum of the elements in the vector MUST be equal to 1.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param depth_choice String. This parameter specifies the choice of the depth.
#' @return A list with the following variables:
#' \item{\code{data}}{A matrix of dimensions #mesh nodes-by-#functions containing the data used in the algorithm.}
#' \item{\code{order}}{Order of the finite elements given as input in IFD.FEM().}
#' \item{\code{weights}}{Weights given as input in IFD.FEM().}
#' \item{\code{ifd}}{IFD computed.}
#' \item{\code{depth}}{Depth computed.}
#' @description This function implements the formula to compute the integrated functional depth of a set of complex multidimensional functional data.
#' The computation relies only on the C++ implementation of the algorithm.
#' @usage IFD.FEM(data, FEMbasis, weights, search = "tree", depth_choice)
#' @export
#' @examples
#' library(fdaPDE)
#' ## example still to be implemented

IFD.FEM <- function(data, FEMbasis, weights, search = "tree", depth_choice) 
{ 
  if(class(FEMbasis$mesh) == "mesh.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }

  # Search algorithm
  if(search=="naive"){
    search=1
  }else if(search=="tree"){
    search=2
  }else if(search=="walking" & class(FEMbasis$mesh) == "mesh.2.5D"){
  stop("walking search is not available for mesh class mesh.2.5D.")
  }else if(search=="walking" & class(FEMbasis$mesh) != "mesh.2.5D"){
    search=3
  }else{
    stop("'search' must must belong to the following list: 'naive', 'tree' or 'walking'.")
  }

  ###################### Checking parameters, sizes and conversion #################################
  checkParametersIFD(data, FEMbasis, weights, search, depth_choice) 
  
  ## Coverting to format for internal usage
  data = as.matrix(data)
  weights = as.vector(weights)
  
  checkParametersSizeIFD(data, FEMbasis, weights) 
  ###################### End checking parameters, sizes and conversion #############################
  
  ###################### C++ Code Execution #########################################################
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){	  
    
    bigsol = CPP_FEM.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)
    
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    
    bigsol = CPP_FEM.manifold.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)
    
  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
    bigsol = CPP_FEM.volume.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)
  }
  
  ###################### Collect Results ############################################################  
  
  data    = bigsol[[1]]
  order   = bigsol[[2]]
  weights = bigsol[[3]]
  ifd     = bigsol[[4]]
  depth   = bigsol[[5]]
  
  reslist = list(data = data, order = order, weights = weights, ifd = ifd, depth = depth)
  return(reslist)
}
