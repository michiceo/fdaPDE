checkParametersIFD <- function(data, FEMbasis, weights, search) 
{
  #################### Parameter Check #########################
  if (is.null(data)) 
    stop("'data' required;  is NULL.")
  else{
    if(any(is.na(data)))
      stop("Missing values not admitted in 'data'.")
  }
  
  if (is.null(FEMbasis)) 
    stop("'FEMbasis' required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'.")
  
  if (is.null(weights))  
    stop("'weights' required;  is NULL.")
  else{
    for(i in 1:length(weights)){
      if(weights[i]<=0)
        stop("'weights' has to have positive members.")
    }
  }
  
  if (is.null(direction_method)) 
    stop("'direction_method' is required;  is NULL.")
  else{
    if(direction_method!="Gradient" && direction_method!="BFGS")
      stop("'direction_method' needs to be either 'Gradient' or 'BFGS'.")
  }

  if(!is.numeric(search))
    stop("'search' needs to be an integer.")
  
  
}


checkParametersSizeIFD <- function(data, FEMbasis, weights) 
{
  if(nrow(data) < 1)
    stop("functions must contain at least one evaluation.")
  if(ncol(data) < 1)                                       
    stop("data must contain at least one function.")
  if(length(weights) != nrow(data))
    stop("The length of weights has to be equal to the number of values")
}