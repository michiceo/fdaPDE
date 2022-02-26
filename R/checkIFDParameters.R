checkParametersIFD <- function(data, FEMbasis, search, depth_choice)
{
  #################### Parameter Check #########################
  if (is.null(data))
    stop("'data' required;  is NULL.")

  if (is.null(FEMbasis))
    stop("'FEMbasis' required;  is NULL.")
  if(class(FEMbasis)!= "FEMbasis")
    stop("'FEMbasis' is not class 'FEMbasis'.")

  if(!is.numeric(search))
    stop("'search' needs to be an integer.")

  if (is.null(depth_choice))
    stop("'depth_choice' is required;  is NULL.")

}


checkParametersSizeIFD <- function(data, FEMbasis, weights)
{
  if(nrow(data) < 1)
    stop("functions must contain at least one evaluation.")
  if(ncol(data) < 1)
    stop("data must contain at least one function.")
}
