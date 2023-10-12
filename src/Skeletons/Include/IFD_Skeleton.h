#ifndef __IFD_SKELETON_H__
#define __IFD_SKELETON_H__

#include "../../FE_Assemblers_Solvers/Include/Finite_Element.h"
#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include "../../Mesh/Include/Mesh.h"
#include "../../FE_Assemblers_Solvers/Include/Matrix_Assembler.h"
#include "../../Global_Utilities/Include/Solver_Definitions.h"

#include "../../Integrated_Functional_Depth/Include/IFD_Integration.h"
#include "../../Integrated_Functional_Depth/Include/IFD_Functional_Boxplot.h"

template<UInt ORDER, UInt mydim, UInt ndim>
SEXP IFD_Skeleton(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string & depth_choice)
{
	// Construct data problem object
	DepthIntegration<ORDER, mydim, ndim> depthIntegration(Rdata, Rorder, Rweights, Rsearch, Rmesh, depth_choice);
	FunctionalBoxplot<ORDER, mydim, ndim> functionalBoxplot(depthIntegration);

	// Perform the whole task
	//functionalBoxplot.apply();
	functionalBoxplot.calculateFunctionalBoxplot();

	// Collect results
	VectorXr ifd = functionalBoxplot.getIFD();
	VectorXr median = functionalBoxplot.getMedian();
	VectorXr firstQuartile = functionalBoxplot.getFirstQuartile();
	VectorXr thirdQuartile = functionalBoxplot.getThirdQuartile();
	VectorXr lowerWhisker = functionalBoxplot.getLowerWhisker();
	VectorXr upperWhisker = functionalBoxplot.getUpperWhisker();
	VectorXr signDepth = functionalBoxplot.getSignDepth();


	// Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 10));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, depthIntegration.dataRows(), depthIntegration.dataCols()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(INTSXP, 1));
	SET_VECTOR_ELT(result, 2, Rf_allocMatrix(REALSXP, depthIntegration.dataRows(), depthIntegration.dataCols()));
	SET_VECTOR_ELT(result, 3, Rf_allocVector(REALSXP, ifd.size()));
	SET_VECTOR_ELT(result, 4, Rf_allocVector(REALSXP, median.size()));
	SET_VECTOR_ELT(result, 5, Rf_allocVector(REALSXP, firstQuartile.size()));
	SET_VECTOR_ELT(result, 6, Rf_allocVector(REALSXP, thirdQuartile.size()));
	SET_VECTOR_ELT(result, 7, Rf_allocVector(REALSXP, lowerWhisker.size()));
	SET_VECTOR_ELT(result, 8, Rf_allocVector(REALSXP, upperWhisker.size()));
	SET_VECTOR_ELT(result, 9, Rf_allocVector(REALSXP, signDepth.size()));



	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < depthIntegration.dataCols(); j++)
	{
		for(UInt i = 0; i < depthIntegration.dataRows(); i++)
			rans[i + depthIntegration.dataRows()*j] = depthIntegration.data()(i, j);
	}

	int *rans1 = INTEGER(VECTOR_ELT(result, 1));
	*rans1 = depthIntegration.getOrder();

	Real *rans2 = REAL(VECTOR_ELT(result, 2));
	for(UInt j = 0; j < depthIntegration.dataCols(); j++)
	{
		for(UInt i = 0; i < depthIntegration.dataRows(); i++)
			rans2[i + depthIntegration.dataRows()*j] = depthIntegration.getWeights()(i, j);
	}

	Real *rans3 = REAL(VECTOR_ELT(result, 3));
	for(UInt i = 0; i < ifd.size(); i++)
	{
		rans3[i] = ifd[i];
	}

	Real *rans4 = REAL(VECTOR_ELT(result, 4));
	for(UInt i = 0; i < median.size(); i++)
	{
		rans4[i] = median[i];
	}

	Real *rans5 = REAL(VECTOR_ELT(result, 5));
	for(UInt i = 0; i < firstQuartile.size(); i++)
	{
		rans5[i] = firstQuartile[i];
	}

	Real *rans6 = REAL(VECTOR_ELT(result, 6));
	for(UInt i = 0; i < thirdQuartile.size(); i++)
	{
		rans6[i] = thirdQuartile[i];
	}

	Real *rans7 = REAL(VECTOR_ELT(result, 7));
	for(UInt i = 0; i < lowerWhisker.size(); i++)
	{
		rans7[i] = lowerWhisker[i];
	}

	Real *rans8 = REAL(VECTOR_ELT(result, 8));
	for(UInt i = 0; i < upperWhisker.size(); i++)
	{
		rans8[i] = upperWhisker[i];
	}

	Real *rans9 = REAL(VECTOR_ELT(result, 9));
	for(UInt i = 0; i < signDepth.size(); i++)
	{
		rans9[i] = signDepth[i];
	}

	UNPROTECT(1);

	return(result);
};

#endif /* __IFD_SKELETON_H__ */
