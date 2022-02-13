#ifndef __IFD_DATA_PROBLEM_H__
#define __IFD_DATA_PROBLEM_H__

#include <string>
#include <utility>
#include <numeric>

#include "../../FdaPDE.h"
#include "IFD_Data.h" // for R/C++ interface
#include "IFD_Depth.h"
#include "IFD_Depth_Factory.h"
#include "../../FE_Assemblers_Solvers/Include/Integration.h"

// This file contains data informations for the computation of the Integrated Functional Depth

template<UInt ORDER, UInt mydim, UInt ndim>
class DataProblem {
private:
	using Integrator = typename FiniteElement<ORDER, mydim, ndim>::Integrator; // using Integrator = typename SpaceIntegratorHelper::Integrator<ORDER, mydim>;
	static constexpr UInt EL_NNODES = how_many_nodes(ORDER, mydim);
	IFDData ifdData_; // for R/C++ interface
	MeshHandler<ORDER, mydim, ndim> mesh_;
	std::shared_ptr<Depth> depth_;
	std::string d_tag; // tag in order to choose the depth
	Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES> PsiQuad_;

	//! A method to compute the matrix which evaluates the basis function at the quadrature Integrator::NNODES
	void fillPsiQuad();

public:
	//! A constructor: it delegates IFDData and MeshHandler constructors.
	DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string& d);

	//! A method to compute the weighted integral of the depth referred to the data.
	const VectorXr FEintegrate_depth(const MatrixXr& X) const;

	//! A method to compute the integral of the weights.
	const VectorXr FEintegral_weights() const;

	// Getters
	//! A method to access the data. It calls the same method of IFDData class.
	inline const MatrixXr & data() const {return ifdData_.data();}
	//! A method to access a specific function. It calls the same method of IFDData class.
	inline const VectorXr & dataC(UInt i) const {return ifdData_.dataC(i);}
	//! A method to access a specific row of the data matrix (evaluations at a specific point). It calls the same method of IFDData class.
	inline const VectorXr & dataR(UInt i) const {return ifdData_.dataR(i);}
	// ! A method to access the number of rows of the data matrix. It calls the same method of IFDData class.
	inline UInt dataRows() const {return ifdData_.dataRows();}
	// ! A method to access the number of columns of the data matrix. It calls the same method of IFDData class.
	inline UInt dataCols() const {return ifdData_.dataCols();}
	//! A method returning the the input order. It calls the same method of IFDData class.
	inline UInt getOrder() const {return ifdData_.getOrder();}
	//! A method returning the weights for the integration. It calls the same method of IFDData class.
	inline const MatrixXr & getWeights() const {return ifdData_.getWeights();} //VectorXr
	//! A method returning the depth of the data (no integration).
	inline const VectorXr getDepth() const {return depth_->depth();}
	inline const VectorXr getDepth(UInt i) const {return depth_->depth(i);}

	// Getters for mesh
	//! A method returning the mesh.
	inline const MeshHandler<ORDER, mydim, ndim>& getMesh() const {return mesh_;}
	// Getter for specific mesh features
	//! A method returning the number of mesh EL_NNODES. It calls the same method of MeshHandler class.
	inline UInt getNumNodes() const {return mesh_.num_nodes();}
	//! A method returning the number of mesh elements. It calls the same method of MeshHandler class.
	inline UInt getNumElements() const {return mesh_.num_elements();}
	//! A method returning a node. It calls the same method of MeshHandler class.
	inline Point<ndim> getPoint(Id id) const {return mesh_.getPoint(id);}
	//! A method returning an element. It calls the same method of MeshHandler class.
	inline Element<EL_NNODES,mydim,ndim> getElement(Id id) const {return mesh_.getElement(id);}
	//! A method returning the element in which the point in input is located. It calls the same method of MeshHandler class.
	inline Element<EL_NNODES,mydim,ndim> findLocation(const Point<ndim>& point) const {return mesh_.findLocation(point);}

	// Getters for matrices
	//! A method returning the PsiQuad_ matrix.
	inline const Eigen::Matrix<Real, Integrator::NNODES, EL_NNODES>& getPsiQuad() const {return PsiQuad_;}
};

#include "IFD_Data_Problem_imp.h"

#endif /*__IFD_DATA_PROBLEM_H__*/
