#ifndef __FE_INTEGRATED_FUNCTIONAL_DEPTH_H__
#define __FE_INTEGRATED_FUNCTIONAL_DEPTH_H__

#include "Data_Problem.h"

//! @brief A class to compute both the classical depth (no integration) and the IFD.

template<UInt ORDER, UInt mydim, UInt ndim>
class FEIFD {
private:
	// A member to acess data problem methods
	const DataProblem<ORDER, mydim, ndim>& dataProblem_;
	VectorXr ifd;
	VectorXr depth;

public:
	//! A constructor
	FEIFD(const DataProblem<ORDER, mydim, ndim>& dP): dataProblem_(dP){};

	//! A method to perform the whole computation of the i.f.d. task.
	void apply();

	// Getters
	//! A method returning the computed i.f.d.
	const VectorXr getIFD() const {return ifd;}
	//! A method returning the classical depth of the data (no integration)
	const VectorXr getDepth() const {return depth;}
};

#include "FE_Integrated_Functional_Depth_imp.h"

#endif /* __FE_INTEGRATED_FUNCTIONAL_DEPTH_H__ */

