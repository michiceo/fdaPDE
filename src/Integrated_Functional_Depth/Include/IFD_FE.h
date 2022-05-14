#ifndef __IFD_FE_H__
#define __IFD_FE_H__

#include "IFD_Data_Problem.h"

//! @brief A class to compute both the classical depth (no integration) and the IFD.

template<UInt ORDER, UInt mydim, UInt ndim>
class FEIFD {
private:
	// A member to acess data problem methods
	const DataProblem<ORDER, mydim, ndim>& dataProblem_;

	VectorXr ifd;
	VectorXr depth;
	VectorXr median;
	VectorXr firstQuartile;
	VectorXr thirdQuartile;
	VectorXr lowerWhisker;
	VectorXr upperWhisker;

public:
	//! A constructor
	FEIFD(const DataProblem<ORDER, mydim, ndim>& dP): dataProblem_(dP){
		ifd = dataProblem_.FEintegrate_depth(dataProblem_.data());
		depth = dataProblem_.getDepth();};

	//! A method to perform the whole computation of the i.f.d. task.
	//void apply();

	void calculateFunctionalBoxplot();

	// Getters
	//! A method returning the computed i.f.d.
	const VectorXr getIFD() const {return ifd;}
	//! A method returning the classical depth of the data (no integration)
	const VectorXr getDepth() const {return depth;}
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
    //! A method returning the data
	const VectorXr getData() const {return dataProblem_.data();}
};

#include "IFD_FE_imp.h"

#endif /* __IFD_FE_H__ */

