#ifndef __IFD_INTEGRATION_IMP_H__
#define __IFD_INTEGRATION_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
DepthIntegration<ORDER, mydim, ndim>::DepthIntegration(SEXP Rdata, SEXP Rorder, SEXP Rweights, SEXP Rsearch, SEXP Rmesh, const std::string& d):
ifdData_(Rdata, Rorder, Rweights), mesh_(Rmesh, INTEGER(Rsearch)[0]), d_tag(d)
{
	fillPsiQuad();
}

template<UInt ORDER, UInt mydim, UInt ndim>
void
DepthIntegration<ORDER, mydim, ndim>::fillPsiQuad()
{
	for(UInt i=0; i<Integrator::NNODES; ++i)
	   PsiQuad_.row(i)=reference_eval_point<EL_NNODES, mydim>(Integrator::NODES[i]);
}

template<UInt ORDER, UInt mydim, UInt ndim>
const VectorXr
DepthIntegration<ORDER, mydim, ndim>::integrate_depth(const MatrixXr& X) const
{
	Real time1 = omp_get_wtime();

	omp_set_num_threads(2);

	using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1>>;

	VectorXr total_sum = VectorXr::Zero(X.cols());

	#pragma omp parallel for default(none) shared(total_sum, X)
	for(UInt triangle = 0; triangle < mesh_.num_elements(); ++triangle){

		MatrixXr X_cap = MatrixXr::Zero(Integrator::NNODES, 1);
		std::unique_ptr<Depth> depth_cap;

		Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);

		MatrixXr sub_w = MatrixXr::Zero(EL_NNODES, X.cols());

		for(Eigen::Index j=0; j < X.cols(); ++j){
			const VectorXr& x = X.col(j);

			Eigen::Matrix<Real, EL_NNODES, 1> sub_x;
			for(UInt i = 0; i < EL_NNODES; ++i){
				sub_x[i] = x[tri_activated[i].getId()];
			  sub_w(i, j) = this->getWeights()(tri_activated[i].getId(), j);
			}

			Eigen::Matrix<Real, Integrator::NNODES, 1> x_cap = (PsiQuad_*sub_x).array();
			X_cap.col(X_cap.cols()-1) = x_cap; // fill the matrix with the transformed functions (due to quadrature formula evaluation)

			if(x != X.col( X.cols()-1 )) // O(1)
				X_cap.conservativeResize(X_cap.rows(), X_cap.cols()+1);
		}

		depth_cap = Depth_factory::createDepth(X_cap, d_tag);

		for(Eigen::Index j=0; j < X.cols(); ++j){
			Eigen::Matrix<Real, Integrator::NNODES, 1> weights = (PsiQuad_*sub_w.col(j)).array();
			Eigen::DiagonalMatrix<Real, Integrator::NNODES, Integrator::NNODES> weighdiag;
			weighdiag.diagonal() = weights;
			#pragma imp atomic update
			total_sum[j] += ( weighdiag * depth_cap->depth(j) ).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
		}
	}
	return total_sum;
}

#endif /* __IFD_INTEGRATION_IMP_H__ */
