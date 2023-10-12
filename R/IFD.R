#' Integrated depth for functions over complicated multidimensional domains
#' @param data A matrix of dimensions #mesh nodes-by-#functions. Data are functions:
#' each row corresponds to the evaluations of the functions at one specific node of the mesh,
#' each column corresponds to the evaluation of a specific function at all the mesh nodes.
#' @param FEMbasis A \code{FEMbasis} object describing the Finite Element basis,
#' as created by \code{\link{create.FEM.basis}}.
#' @param weights A weight function. The integral of the function MUST be equal to 1.
#' @param search a flag to decide the search algorithm type (tree or naive or walking search algorithm).
#' @param depth_choice String. This parameter specifies the choice of the depth.
#' @param plot Boolean. This parameter specifies whether one wants the plot of the functional boxplot or not.
#' @param func Vector. This parameter specifies the function of the dataset whose location in the functional boxplot will be shown.
#' @return A list with the following variables:
#' \item{\code{data}}{A matrix of dimensions #mesh nodes-by-#functions containing the data used in the algorithm.}
#' \item{\code{order}}{Order of the finite elements given as input in IFD.FEM().}
#' \item{\code{weights}}{Evaluation of weight function at nodes.}
#' \item{\code{ifd}}{Weighted depth (Integrated Functional Depth).}
#' \item{\code{median, firstQuartile, thirdQuartile, lowerWhisker, upperWhisker}}{The quantities to visualize the "functional boxplot".}
#' @description This function implements the formula to compute the integrated depth for a set of functions over complicated multidimensional domains, with the option to return also the "functional boxplot".
#' The computation relies only on the C++ implementation of the algorithm.
#' @usage IFD.FEM(data, FEMbasis, weights = NULL, search = "tree", depth_choice, plot = TRUE, func = data[,1])
#' @export
#' @examples
#' library(fdaPDE)
#'
#' ## Create a 2D mesh over a squared domain
#' x = seq(0,1, length.out = 3)
#' y = x
#' locations = expand.grid(x,y)
#' mesh = create.mesh.2D(locations)
#' plot(mesh)
#' nnodes = dim(mesh$nodes)[1]
#' FEMbasis = create.FEM.basis(mesh)
#'
#' ## Generate data
#' data = NULL
#' for(ii in 1:50){
#'   a1 = rnorm(1, mean = 1, sd = 1)
#'   a2 = rnorm(1, mean = 1, sd = 1)
#'
#'   func_evaluation = numeric(nrow(mesh$nodes))
#'   for (i in 0:(nrow(mesh$nodes)-1)){
#'     func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +
#'                            a2* sin(2*pi*mesh$nodes[i+1,2])+ 1
#'   }
#'   datum = func_evaluation + rnorm(nrow(mesh$nodes), mean = 0, sd = 0.5)
#'   data = cbind(data, datum)
#'   colnames(data) = NULL
#' }
#'
#' ## Computation of the depth
#' sol <- IFD.FEM(data = data, FEMbasis = FEMbasis, depth_choice = "MHRD")
#'

IFD.FEM <- function(data, FEMbasis, weights = NULL, search = "tree", depth_choice, plot = FALSE, func = 0)
{
  if(class(FEMbasis$mesh) == "mesh.1.5D"){
    ndim = 2
    mydim = 1 
  }else if(class(FEMbasis$mesh) == "mesh.2D"){
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
  }else if(search=="walking" && class(FEMbasis$mesh) == "mesh.2.5D"){
    stop("walking search is not available for mesh class mesh.2.5D.")
  }else if(search=="walking" && class(FEMbasis$mesh) != "mesh.2.5D"){
    search=3
  }else{
    stop("'search' must belong to the following list: 'naive', 'tree' or 'walking'.")
  }

  ###################### Checking parameters, sizes and conversion #################################

  if(is.null(weights)){
    n <- dim(data)[2]
    p <- dim(data)[1]

    w <- function(nfun, npoints){
      phi_num <- nfun - rowSums(is.na(data))
      phi_den  <- apply(data, 2, function(x) sum(phi_num[!is.na(x)]))

      output <- matrix(0, npoints, nfun)
      for(i in 1:nfun){
        for(j in 1:npoints){
          if(is.na(data[j,i])){
            output[j,i] <- 0
          }
          else{
            output[j,i] <- (phi_num[j]*npoints)/(phi_den[i])
          }
        }
      }

      output # matrix with the weights for each function (considers the NA)
    }

    weights <- w(n, p)
  }
  else{
    w <- weights(FEMbasis$mesh$nodes)
  }

  
  checkParametersIFD(data, FEMbasis, search, depth_choice)

  ## Coverting to format for internal usage
  data = as.matrix(data)
  weights = as.matrix(weights)

  checkParametersSizeIFD(data, FEMbasis)
  ###################### End checking parameters, sizes and conversion #############################

  ###################### C++ Code Execution #########################################################
  bigsol = NULL
  if(class(FEMbasis$mesh) == 'mesh.1.5D'){

    bigsol = CPP_FEM_1.5D.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)

  } else if(class(FEMbasis$mesh) == 'mesh.2D'){

    bigsol = CPP_FEM.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)

  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){

    bigsol = CPP_FEM.manifold.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)

  } else if(class(FEMbasis$mesh) == 'mesh.3D'){
   
    bigsol = CPP_FEM.volume.IFD(data, FEMbasis, ndim, mydim, weights, search, depth_choice)
  }

  ###################### Collect Results ############################################################

  data     = bigsol[[1]]
  order    = bigsol[[2]]
  weights  = bigsol[[3]]
  ifd      = bigsol[[4]]
  median   = bigsol[[5]]
  firstQuartile   = bigsol[[6]]
  thirdQuartile   = bigsol[[7]]
  lowerWhisker    = bigsol[[8]]
  upperWhisker    = bigsol[[9]]
  signDepth       = bigsol[[10]]

  if(plot){
    if(class(FEMbasis$mesh) == "mesh.2D"){
      max<-max(data)
      min<-min(data)
      par(mfrow=c(1, 5))
      plot.image.2D(FEM(lowerWhisker, FEMbasis), max, min)
      plot.image.2D(FEM(firstQuartile, FEMbasis), max, min)
      plot.image.2D(FEM(median, FEMbasis), max, min)
      plot.image.2D(FEM(thirdQuartile, FEMbasis), max, min)
      plot.image.2D(FEM(upperWhisker, FEMbasis), max, min)

      par(mfrow=c(1, 5))
      plot.image.diff_lw_q1.2D(FEM(func, FEMbasis), FEM(lowerWhisker, FEMbasis ), max, min)
      plot.image.diff_lw_q1.2D(FEM(func, FEMbasis), FEM(firstQuartile, FEMbasis ), max, min)
      plot.image.2D(FEM(func, FEMbasis), max, min)
      plot.image.diff_uw_q3.2D(FEM(func, FEMbasis), FEM(thirdQuartile, FEMbasis ), max, min)
      plot.image.diff_uw_q3.2D(FEM(func, FEMbasis), FEM(upperWhisker, FEMbasis ), max, min)
    }else{
      max<-max(data[!is.na(data)]) 
      min<-min(data[!is.na(data)])

      plot.FEM(FEM(median, FEMbasis), max, min)
      plot.FEM(FEM(firstQuartile, FEMbasis), max, min)
      plot.FEM(FEM(thirdQuartile, FEMbasis), max, min)
      plot.FEM(FEM(lowerWhisker, FEMbasis), max, min)
      plot.FEM(FEM(upperWhisker, FEMbasis), max, min)

      if(func[1] != 0){
        plot.diff_lw_q1.3D(FEM(func, FEMbasis), FEM(firstQuartile, FEMbasis), max, min)
        plot.diff_uw_q3.3D(FEM(func, FEMbasis), FEM(thirdQuartile, FEMbasis), max, min)
        plot.diff_lw_q1.3D(FEM(func, FEMbasis), FEM(lowerWhisker, FEMbasis), max, min)
        plot.diff_uw_q3.3D(FEM(func, FEMbasis), FEM(upperWhisker, FEMbasis), max, min)
        plot(FEM(func, FEMbasis))
      }
    }
  }

  # if (ifd == 0. && depth==0.)
    # stop("integral of weight function != 1. Give another weight function.")

  reslist = list(data = data, order = order, weights = weights, ifd = ifd, median = median,
    firstQuartile = firstQuartile, thirdQuartile = thirdQuartile, lowerWhisker = lowerWhisker, upperWhisker = upperWhisker, signDepth = signDepth)
  return(reslist)
}
