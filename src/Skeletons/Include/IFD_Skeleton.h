#ifndef __IFD_SKELETON_H__
#define __IFD_SKELETON_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

#include "../../Integrated_Functional_Depth/Include/IFD_Data_Problem.h"
#include "../../Integrated_Functional_Depth/Include/IFD_FE.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP IFD_Skeleton(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string & depth_choice)
{
	// Construct data problem object
	DataProblem<ORDER, mydim, ndim> dataProblem(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice);
	FEIFD<ORDER, mydim, ndim> feifd(dataProblem);

	// Perform the whole task
	//feifd.apply();
	feifd.calculateFunctionalBoxplot();

	// Collect results
	VectorXr ifd = feifd.getIFD();
	VectorXr depth = feifd.getDepth();
	VectorXr median = feifd.getMedian();
	VectorXr firstQuartile = feifd.getFirstQuartile();
	VectorXr thirdQuartile = feifd.getThirdQuartile();
	VectorXr lowerWhisker = feifd.getLowerWhisker();
	VectorXr upperWhisker = feifd.getUpperWhisker();

	// Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 10));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, dataProblem.dataRows(), dataProblem.dataCols()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(INTSXP, 1));
	SET_VECTOR_ELT(result, 2, Rf_allocMatrix(REALSXP, dataProblem.dataRows(), dataProblem.dataCols()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, ifd.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, depth.size()));
	SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, median.size()));
	SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, firstQuartile.size()));
	SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, thirdQuartile.size()));
	SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, lowerWhisker.size()));
	SET_VECTOR_ELT(result, 9, Rf_allocVector(REALSXP, upperWhisker.size()));


	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < dataProblem.dataCols(); j++)
	{
		for(UInt i = 0; i < dataProblem.dataRows(); i++)
			rans[i + dataProblem.dataRows()*j] = dataProblem.data()(i, j);
	}

	int *rans1 = INTEGER(VECTOR_ELT(result, 1));
	*rans1 = dataProblem.getOrder();

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt j = 0; j < dataProblem.dataCols(); j++)
	{
		for(UInt i = 0; i < dataProblem.dataRows(); i++)
			rans2[i + dataProblem.dataRows()*j] = dataProblem.getWeights()(i, j);
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

	Real *rans5 = REAL(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < median.size(); i++)
	{
		rans5[i] = median[i];
	}

	Real *rans6 = REAL(VECTOR_ELT(result, 6));
	for(UInt i = 0; i < firstQuartile.size(); i++)
	{
		rans6[i] = firstQuartile[i];
	}

	Real *rans7 = REAL(VECTOR_ELT(result, 7));
	for(UInt i = 0; i < thirdQuartile.size(); i++)
	{
		rans7[i] = thirdQuartile[i];
	}

	Real *rans8 = REAL(VECTOR_ELT(result, 8));
	for(UInt i = 0; i < lowerWhisker.size(); i++)
	{
		rans8[i] = lowerWhisker[i];
	}

	Real *rans9 = REAL(VECTOR_ELT(result, 9));
	for(UInt i = 0; i < upperWhisker.size(); i++)
	{
		rans9[i] = upperWhisker[i];
	}

	UNPROTECT(1);

	return(result);
};

#endif /* __IFD_SKELETON_H__ */
