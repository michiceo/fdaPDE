#ifndef __IFD_FUNCTIONAL_BOXPLOT_H__
#define __IFD_FUNCTIONAL_BOXPLOT_H__

#include "IFD_Integration.h"

//! @brief A class to compute both the classical depth (no integration) and the IFD.

template<UInt ORDER, UInt mydim, UInt ndim>
class FunctionalBoxplot {
private:
	// A member to acess data problem methods
	const DepthIntegration<ORDER, mydim, ndim>& depthIntegration_;

	VectorXr ifd;
	VectorXr median;
	VectorXr firstQuartile;
	VectorXr thirdQuartile;
	VectorXr lowerWhisker;
	VectorXr upperWhisker;
	VectorXr signDepth;

public:
	//! A constructor
	FunctionalBoxplot(const DepthIntegration<ORDER, mydim, ndim>& dI): depthIntegration_(dI){
		ifd = depthIntegration_.integrate_depth(depthIntegration_.data());
	}

	void calculateFunctionalBoxplot();

	// Getters
	//! A method returning the computed i.f.d.
	const VectorXr getIFD() const {return ifd;}
	//! A method returning the median
	const VectorXr getMedian() const {return median;}
	//! A method returning the first quartile
	const VectorXr getFirstQuartile() const {return firstQuartile;}
	//! A method returning the third quartile
	const VectorXr getThirdQuartile() const {return thirdQuartile;}
	//! A method returning the lower whisker
	const VectorXr getLowerWhisker() const {return lowerWhisker;}
	//! A method returning the upper whisker
	const VectorXr getUpperWhisker() const {return upperWhisker;}
	//! A method returning the signDepth
	const VectorXr getSignDepth() const {return signDepth;}
    //! A method returning the data
	const VectorXr getData() const {return depthIntegration_.data();}
};

#include "IFD_Functional_Boxplot_imp.h"

#endif /* __IFD_FUNCTIONAL_BOXPLOT_H__ */
