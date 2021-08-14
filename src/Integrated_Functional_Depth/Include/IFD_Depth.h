#ifndef __IFD_DEPTH_H__
#define __IFD_DEPTH_H__

#include <algorithm>
#include <string>
#include "../../FdaPDE.h"

// This file contains the possible functional depth measures choices
// for the computation of the Integrated Functional Depth

/*! @brief An abstract base class dealing with the depth computation.
*/

class Depth{

public:
	//! A Constructor.
	Depth(const MatrixXr& m): m_(m){};
	//! A Destructor.
	virtual ~Depth(){};
	//! A pure virtual method to compute the depth of all data.
	virtual const VectorXr compute_depth() const = 0;
	//! A pure virtual method to compute the depth of the k-th function at each mesh node.
	virtual const VectorXr compute_depth(UInt k) const = 0;

protected:
	//! A matrix of data
	MatrixXr m_;
	// number of functions
	const UInt n_ = m_.cols();
	// number of points
	const UInt p_ = m_.rows();

};

/*! @brief A class dealing with the computation of the Modified Half Region Depth.
*/
class MHRD: public Depth{

public:
	//! A Constructor
	MHRD(const MatrixXr& m);
	//! An overridden method to compute the depth chosen of all data.
	const VectorXr compute_depth() const override;
	//! An overridden method to compute the depth chosen of the k-th function at each mesh node.
	const VectorXr compute_depth(UInt k) const override;

private:
	//! A method to perform the ranking of the elements in the vector.
	const VectorXi ranking(const VectorXr& v) const;

};

#endif /* __IFD_DEPTH_H__ */

