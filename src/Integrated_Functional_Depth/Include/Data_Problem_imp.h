#ifndef __DATA_PROBLEM_IMP_H__
#define __DATA_PROBLEM_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
DataProblem<ORDER, mydim, ndim>::DataProblem(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string& d):
ifdData_(Rdata, Rorder, Rweights), mesh_(Rmesh, INTEGER(Rsearch)[0]), d_tag(d)
{
	// depth_ = Depth_factory::createDepth(this->data(), d_tag);
	fillPsiQuad();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void
DataProblem<ORDER, mydim, ndim>::fillPsiQuad()
{
	for(UInt i=0; i<Integrator::NNODES; ++i)
	   PsiQuad_.row(i)=reference_eval_point<EL_NNODES, mydim>(Integrator::NODES[i]);
}

template<UInt ORDER, UInt mydim, UInt ndim>
const VectorXr
DataProblem<ORDER, mydim, ndim>::FEintegrate_depth(const MatrixXr& X) const
{
	using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1>>;

	VectorXr total_sum = VectorXr::Zero(X.cols());

	for(UInt triangle = 0; triangle < mesh_.num_elements(); ++triangle){

		MatrixXr X_cap = MatrixXr::Zero(Integrator::NNODES, 1);
		std::shared_ptr<Depth> depth_cap;

		Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);
		
		Eigen::Matrix<Real, EL_NNODES, 1> sub_w;
    		for(UInt i = 0; i < EL_NNODES; ++i){
			sub_w[i] = this->getWeights()[tri_activated[i].getId()];
		}

		for(Eigen::Index j=0; j < X.cols(); ++j){
			const VectorXr& x = X.col(j);

			Eigen::Matrix<Real, EL_NNODES, 1> sub_x;
			for(UInt i = 0; i < EL_NNODES; ++i){
				sub_x[i] = x[tri_activated[i].getId()];
			}

			Eigen::Matrix<Real, Integrator::NNODES, 1> x_cap = (PsiQuad_*sub_x).array();
			X_cap.col(X_cap.cols()-1) = x_cap; // fill the matrix with the transformed functions (due to quadrature formula evaluation)
			
			if(x != X.col( X.cols()-1 ))
				X_cap.conservativeResize(X_cap.rows(), X_cap.cols()+1);
		}

		depth_cap = Depth_factory::createDepth(X_cap, d_tag);
		Eigen::Matrix<Real, Integrator::NNODES, 1> weights = (PsiQuad_*sub_w).array();

		for(Eigen::Index i = 0; i<X.cols(); ++i)
			total_sum[i] += depth_cap->compute_depth()[i] * weights.dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
	}

	return total_sum;
}

#endif /* __DATA_PROBLEM_IMP_H__ */

