#ifndef __IFD_SKELETON_H__
#define __IFD_SKELETON_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

#include "../../Integrated_Functional_Depth/Include/Data_Problem.h"
#include "../../Integrated_Functional_Depth/Include/FE_Integrated_Functional_Depth.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP IFD_Skeleton(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string & depth_choice)
{
	// Construct data problem object
	DataProblem<ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice);
	
	// Construct FEIFD object
	FEIFD<ORDER, mydim, ndim> feifd(dataProblem);

	// Perform the whole task
	feifd.apply();
	
	// Collect results
	VectorXr ifd = feifd.getIFD();
	VectorXr depth = feifd.getDepth();
	
	// Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 5));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, dataProblem.dataRows(), dataProblem.dataCols()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(INTSXP, 1));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, dataProblem.getWeights().size()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, ifd.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, depth.size()));

	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < dataProblem.dataCols(); j++)
	{
		for(UInt i = 0; i < dataProblem.dataRows(); i++)
			rans[i + dataProblem.dataRows()*j] = dataProblem.data()(i, j);
	}

	int *rans1 = INTEGER(VECTOR_ELT(result, 1));
	*rans1 = dataProblem.getOrder();

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < dataProblem.getWeights().size(); i++)
	{
		rans2[i] = dataProblem.getWeights()[i];
	}
	
	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < ifd.size(); i++)
	{
		rans3[i] = ifd[i];
	}
	
	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < depth.size(); i++)
	{
		rans4[i] = depth[i];
	}

	UNPROTECT(1);

	return(result);
};

#endif /* __IFD_SKELETON_H__ */

