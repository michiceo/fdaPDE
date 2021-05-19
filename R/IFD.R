# still missing R documentation
#' @param 
#' @return 
#' @description 
#' @usage 
#' @export
#' @examples

#(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rweights, SEXP Rsearch)
IFD.FEM <- function(data, FEMbasis, weights, search = "tree") 
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
  checkParametersIFD(data, FEMbasis, weights, search) 
  
  ## Coverting to format for internal usage
  data = as.matrix(data)
  weights = as.vector(weights)
  
  
  checkParametersSizeIFD(data, FEMbasis, weights) 
  ###################### End checking parameters, sizes and conversion #############################
  
  
  ###################### C++ Code Execution #########################################################
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.2D'){	  
    
    bigsol = CPP_FEM.IFD(data, FEMbasis, ndim, mydim, weights, search)
    
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    
    bigsol = CPP_FEM.manifold.IFD(data, FEMbasis, ndim, mydim, weights, search)
    
  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
    bigsol = CPP_FEM.volume.IFD(data, FEMbasis, ndim, mydim, weights, search)
  }
  
  ###################### Collect Results ############################################################  
  
  data    = bigsol[[1]]
  order   = bigsol[[2]]
  weights = bigsol[[3]]
  
  reslist = list(data = data, order = order, weights = weights)
  return(reslist)
}
