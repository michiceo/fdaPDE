#include "../../FdaPDE.h"
#include "../../Skeletons/Include/IFD_Skeleton.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

#include "../Include/IFD_Data_Problem.h"

extern "C" {

    //! This function manages the various options for IFD-PDE algorithm
    /*!
    	This function is then called from R code.
    	\param Rdata an R-matrix containing the data.
    	\param Rmesh an R-object containg the output mesh from Trilibrary
    	\param Rorder an R-integer containing the order of the approximating basis.
    	\param Rmydim an R-integer containing the dimension of the problem we are considering.
    	\param Rndim an R-integer containing the dimension of the space in which the location are.
    	\param Rweights an R-vector containing the weights for the integration.
    	\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
    	\param Rdepth an R-string containing the depth choice.

    	\return R-list containg solutions.
    */

    SEXP Integrated_Functional_Depth(SEXP Rdata, SEXP Rmesh, SEXP Rorder, SEXP Rmydim, SEXP Rndim, SEXP Rweights, SEXP Rsearch, SEXP Rdepth)
    {
        UInt order = INTEGER(Rorder)[0];
        UInt mydim = INTEGER(Rmydim)[0];
        UInt ndim  = INTEGER(Rndim)[0];
        
        std::string depth_choice = CHAR(STRING_ELT(Rdepth, 0));

	if(order== 1 && mydim==2 && ndim==2)
		return(IFD_Skeleton<1, 2, 2>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	else if(order== 2 && mydim==2 && ndim==2)
		return(IFD_Skeleton<2, 2, 2>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	else if(order== 1 && mydim==2 && ndim==3)
		return(IFD_Skeleton<1, 2, 3>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	else if(order== 2 && mydim==2 && ndim==3)
		return(IFD_Skeleton<2, 2, 3>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	else if(order == 1 && mydim==3 && ndim==3)
		return(IFD_Skeleton<1, 3, 3>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	else if(order == 2 && mydim==3 && ndim==3)
		return(IFD_Skeleton<2, 3, 3>(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice));
	
	return(NILSXP);
    }

}
