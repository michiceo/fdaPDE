#ifndef __IFD_SKELETON_H__
#define __IFD_SKELETON_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

#include "../../Integrated_Functional_Depth/Include/IFD_Data.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP IFD_Skeleton(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh)
{
	// Construct data object
	IFDData data(Rdata, Rorder, Rweights);
	
	// Construct mesh object
	MeshHandler<ORDER, mydim, ndim> mesh(Rmesh, INTEGER(Rsearch)[0]);

	// Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 3));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, data.dataRows(), data.dataCols()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(INTSXP, 1));
	SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, data.getWeights().size()));

	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < data.dataCols(); j++)
	{
		for(UInt i = 0; i < data.dataRows(); i++)
			rans[i + data.dataRows()*j] = data.data()(i, j);
	}

	int *rans1 = INTEGER(VECTOR_ELT(result, 1));
	*rans1 = data.getOrder();

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt i = 0; i < data.getWeights().size(); i++)
	{
		rans2[i] = data.getWeights()[i];
	}

	UNPROTECT(1);

	return(result);
};

#endif /* __IFD_SKELETON_H__ */

