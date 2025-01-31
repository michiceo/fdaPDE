#ifndef __DE_DATA_H__
#define __DE_DATA_H__

#include "../../FdaPDE.h"
#include "../../Mesh/Include/Mesh_Objects.h"
#include <array>

// This file contains the R/C++ data conversion for the Density Estimation problem

/*! @brief An IO handler class for objects passed from R.
 * This class, given the spatial data from R, converts them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
template <UInt ndim>
class  DEData{
	private:
		// Data are spatial locations.
		std::vector<Point<ndim>> data_;
		// Finite element order.
		UInt order_;
		// Initial coefficients for the density.
		VectorXr fvec_;
		// Time step parameter for the heat diffusion process.
		Real heatStep_;
		// Number of iterations for the heat diffusion process.
		UInt heatIter_;
		// Penalization parameters (in space). The best one is chosen via k fold cross validation.
		std::vector<Real> lambda_;
		// Number of folds for cross validation.
		UInt Nfolds_;
		// Number of simulations of the optimization algorithm.
		UInt nsim_;
		// Optimization parameters for fixed step methods.
		std::vector<Real> stepProposals_;
		// Tolerances for optimization algorithm termination criteria.
		Real tol1_;
		Real tol2_;
		// A boolean: true if the user wants to see the value of the functional printed during the optimization descent.
		bool print_;
		// Integer specifying the search algorithm type (tree or naive search algorithm),
		UInt search_;

		// Auxiliary methods used in the constructor.
		void setData(SEXP Rdata);
		void setFvec(SEXP Rfvec);
		void setLambda(SEXP Rlambda);
		void setStepProposals(SEXP RstepProposals);

	public:
		// Constructors
		DEData(){};

		explicit DEData(const std::vector<Point<ndim>>& data, const UInt& order, const VectorXr& fvec, Real heatStep,
                        UInt heatIter, const std::vector<Real>& lambda, const UInt& nfolds, const UInt& nsim,
                        const std::vector<Real>& stepProposals, Real tol1, Real tol2, bool print, UInt search);

		/*! Constructor useful for the R/C++ interface.
			It initializes the object storing the R given objects.
			\param Rdata an R-matrix containing the data.
			\param Rorder an R-integer containing the order of the approximating basis.
			\param Rfvec an R-vector containing the the initial solution coefficients given by the user.
			\param RheatStep an R-double containing the step for the heat equation initialization.
			\para, RheatIter an R-integer containing the number of iterations to perfrom the heat equation initialization.
			\param Rlambda an R-vector containing the penalization terms (in space).
			\param Rnfolds an R-integer specifying the number of folds for cross validation.
			\param Rnsim an R-integer specifying the number of iterations to use in the optimization algorithm.
			\param RstepProposals an R-vector containing the step parameters useful for the descent algotihm.
			\param Rtol1 an R-double specifying the tolerance to use for the termination criterion based on the percentage differences.
			\param Rtol2 an R-double specifying the tolerance to use for the termination criterion based on the norm of the gradient.
			\param Rprint and R-integer specifying if print on console.
			\param Rsearch an R-integer to decide the search algorithm type (tree or naive search algorithm).
		*/
		explicit DEData(SEXP Rdata, SEXP Rorder, SEXP Rfvec, SEXP RheatStep, SEXP RheatIter, SEXP Rlambda, SEXP Rnfolds,
                        SEXP Rnsim, SEXP RstepProposals, SEXP Rtol1, SEXP Rtol2, SEXP Rprint, SEXP Rsearch);

		// Getters
		//! A method to access the data.
		std::vector<Point<ndim>>& data() {return data_;}
		//! A const method to access the data.
		const std::vector<Point<ndim>>& data() const {return data_;}
		//! A method to access a datum.
		Point<ndim>& data(UInt i) {return data_[i];}
		//! A const method to access a datum.
		const Point<ndim>& data(UInt i) const {return data_[i];}
		//! A method returning the number of observations.
		UInt dataSize() const {return data_.size();}
		//! A method returning the the input order.
		UInt getOrder() const {return order_;}
		//! A method returning the initial coefficients for the density.
		VectorXr getFvec() const {return fvec_;}
		//! A method returning the heat diffusion process alpha parameter.
		Real getHeatStep() const {return heatStep_;}
		//! A method returning the number of iterations for the heat diffusion process.
		UInt getHeatIter() const {return heatIter_;}
		//! A method returning a bool which says if there is a user's initial density.
		bool isFvecEmpty() const {return fvec_.size() == 0;}
		//! A method returning the penalization parameters (in space).
		Real getLambda(UInt i) const {return lambda_[i];}
		//! A method returning the number of lambdas (in space).
		UInt getNlambda() const {return lambda_.size();}
		//! A method returning the number of folds for CV.
		UInt getNfolds() const {return Nfolds_;}
		//! A method returning the number of iterations to use in the optimization algorithm.
		UInt getNsimulations() const {return nsim_;}
		//! A method returning the number of parameters for fixed step methods.
		UInt getNstepProposals() const {return stepProposals_.size();}
		//! A method returning a parameter for fixed step methods.
		Real getStepProposals(UInt i) const {return stepProposals_[i];}
		//! A method returning the tolerance for optimization algorithm first termination criterion.
		Real getTol1() const {return tol1_;}
		//! A method returning the tolerance for optimization algorithm second termination criterion.
		Real getTol2() const {return tol2_;}
		//! A method returning the boolean print member.
		bool Print() const {return print_;}
		//! A method returning the integer that specifies the search algorithm type.
		UInt getSearch() const {return search_;}

		// Print
		//! A method printing data.
		void printData(std::ostream & out) const;

};


/*! @brief An IO handler class for objects passed from R.
 * This class, given the temporal data from R, converts them in a C++ format, offering a
 * series of method for their access, so isolating the more the possible the specific
 * code for R/C++ data conversion.
*/
// This class is used to manage temporal data in spatio-temporal density estimation problems.
class DEData_time {
private:
    // Time observations during which spatial locations are observed (ordered as they appear in the clean dataset,
    // possibly with duplicates).
    std::vector<Real> data_time_;
    // Time observations in chronological order and without duplicates.
    std::vector<Real> times_;
    // Penalization parameters (in time). The best one is chosen with k fold cross validation.
    std::vector<Real> lambda_time_;
    // Data structure connecting time indices (with respect to the chronological order) to the positions of locations
    // (ordered as data appear in the clean dataset) that are observed at the time instants indexed as above.
    // This structure is useful to build the Upsilon_ FEmatrix in the DataProblem_time class.
    std::vector<std::vector<UInt>> Times2Locations_;

    // Auxiliary methods used in the constructor.
    void setDataTime(SEXP Rdata_time);
    void setLambdaTime(SEXP Rlambda_time);

public:
    // Constructors
    DEData_time(){};

    explicit DEData_time(const std::vector<Real>& data_time, const std::vector<Real>& lambda_time);

    /*! Constructor useful for the R/C++ interface.
        It initializes the object storing the R given objects.
        \param Rdata_time an R-vector containing the time observations.
        \param Rlambda_time an R-vector containing the penalization terms (in time).
    */
    explicit DEData_time(SEXP Rdata_time, SEXP Rlambda_time);

    //! A method to generate the Times2Locations_ and times_ data structures.
    // This method is public because it is called after the data cleaning stage, which is performed in the
    // Data_Problem_time class (after spatial mesh and temporal mesh are provided).
    void setTimes2Locations();

    // Getters
    //! A method to access the data.
    std::vector<Real>& data() {return data_time_;}
    //! A const method to access the data.
    const std::vector<Real>& data() const {return data_time_;}
    //! A const method to access the Times2Locations_ data structure at the i-th time index.
    const std::vector<UInt>& getTimes2Locations(UInt i) const {return Times2Locations_[i];}
    //! A method to access the data (without duplicates).
    std::vector<Real>& times() {return times_;}
    //! A const method to access the data (without duplicates).
    const std::vector<Real>& times() const {return times_;}
    //! A const method to access the i-th time data.
    const Real time(UInt i) const {return times_.size()!=0 ? times_[i] : data_time_[i];}
    //! A method returning the number of distinct time data.
    UInt getNTimes() const {return times_.size()!=0 ? times_.size() : data_time_.size();}
    //! A method returning the number of observations.
    UInt dataSize() const {return data_time_.size();}
    //! A method returning the penalization parameters (in time).
    Real getLambda_time(UInt i) const {return lambda_time_[i];}
    //! A method returning the number of lambdas (in time).
    UInt getNlambda_time() const {return lambda_time_.size();}

    // Print
    //! A method printing data.
    void printData(std::ostream& out) const;
    //! A method printing the Times2Locations_ data structure.
    void printTimes2Locations(std::ostream& out) const;
};

#include "DE_Data_imp.h"

#endif
