#ifndef __IFD_DEPTH_H__
#define __IFD_DEPTH_H__

#include <algorithm>
#include <string>
#include <tuple>
#include <type_traits>
#include "../../FdaPDE.h"

// This file contains the possible functional depth measures choices
// for the computation of the Integrated Functional Depth

/*! @brief An abstract base class dealing with the depth computation.
*/

namespace DepthOutputTypeHelper{
	template<typename Tuple>
	using output_type = typename std::conditional<std::is_same<VectorXr,Tuple>::value,VectorXr,Tuple>::type;
};

class Depth{

public:
	//! A Constructor.
	Depth(const MatrixXr& m): m_(m){};
	//! A Destructor.
	virtual ~Depth(){};
	//! A pure virtual method to compute the depth of the k-th function at each mesh node.
	virtual const VectorXr depth(UInt k) const = 0;

protected:
	//! A matrix of data
	MatrixXr m_;
	// number of functions
	const UInt n_ = m_.cols();
	// number of points
	const UInt p_ = m_.rows();
	//! A method to compute the number of nan in a vector.
	const UInt isnan_vector(const VectorXr& v) const;
	//! A method to perform the ranking of the elements in the vector.
	const VectorXi ranking(const VectorXr& v) const;

};

/*! @brief A class dealing with the computation of the Modified Half Region Depth.
*/
class MHRD: public Depth{

	using output_type = typename DepthOutputTypeHelper::output_type<std::tuple<VectorXr, VectorXr, VectorXr>>;

private:
	//! A method to compute the depth chosen of the k-th function at each mesh node.
	const output_type compute_depth(UInt k) const;

public:
	//! A Constructor
	MHRD(const MatrixXr& m);
	//! Overriden method
	const VectorXr depth(UInt k) const override {return std::get<2>(this->compute_depth(k));};

};

/*! @brief A class dealing with the computation of the Modified Half Region Depth.
*/
class MBD: public Depth{

	using output_type = typename DepthOutputTypeHelper::output_type<VectorXr>;

private:
	//! A method to compute the depth chosen of the k-th function at each mesh node.
	const output_type compute_depth(UInt k) const;

public:
	//! A Constructor
	MBD(const MatrixXr& m);
	//! Overriden method
	const VectorXr depth(UInt k) const override {return this->compute_depth(k);};

};

#endif /* __IFD_DEPTH_H__ */
