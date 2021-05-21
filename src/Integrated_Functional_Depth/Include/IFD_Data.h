#ifndef __IFD_DATA_H__
#define __IFD_DATA_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"

// This file contains the R/C++ data conversion for the Integrated
// Functional Depth (IFD) computation problem

/*! @brief An IO handler class for objects passed from R.
 * This class, given the data from R, convert them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/

class IFDData {
private:
	// Data = evaluation of the functions at each point
	MatrixXr data_;
	// Finte element order
	UInt order_;
	// Weights for the integration of the depth
	VectorXr weights_;

	// Auxiliary methods used in the constructor
	void setData(SEXP Rdata);
	void setWeights(SEXP Rweights);

public:
	// Constructors
	IFDData(){};

	explicit IFDData(const MatrixXr & data, const UInt & order, const VectorXr & weights);

	/*! Costructor useful for the R C++ interface.
				It initializes the object storing the R given objects.
				\param Rdata an R-matrix containing the data.
				\param Rorder an R-integer containing the order of the approximating basis.
				\param Rweights an R-vector containing the weights for the integration.
	*/
	explicit IFDData(SEXP Rdata, SEXP Rorder, SEXP Rweights);

	// Getters
	//! A method to access the data.
	MatrixXr & data() {return data_;}
	//! A const method to access the data.
	const MatrixXr & data() const {return data_;}
	//! A method to access a specific function.
	VectorXr dataC(UInt i) {return data_.col(i);}
	//! A const method to access a specific function.
	const VectorXr dataC(UInt i) const {return data_.col(i);}
	//! A method to access a specific row of the data matrix (evaluations at a specific point).
	VectorXr dataR(UInt i) {return data_.row(i);}
	//! A const method to access a specific row of the data matrix (evaluations at a specific point).
	const VectorXr dataR(UInt i) const {return data_.row(i);}
	// ! A method to access the number of rows of the data matrix.
	UInt dataRows() const {return data_.rows();}
	// ! A method to access the number of columns of the data matrix.
	UInt dataCols() const {return data_.cols();}
	//! A method returning the the input order.
	UInt getOrder() const {return order_;}
	//! A method returning the weights for the integration.
	const VectorXr & getWeights() const {return weights_;}

	// Print
	//! A method printing data.
	void printData(std::ostream & out) const;
};

#endif /* __IFD_DATA_H__ */

