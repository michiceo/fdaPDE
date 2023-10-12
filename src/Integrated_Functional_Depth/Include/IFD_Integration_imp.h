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
	using EigenMap2WEIGHTS = Eigen::Map<const Eigen::Matrix<Real, Integrator::NNODES, 1>>;

	VectorXr total_sum = VectorXr::Zero(X.cols());
	VectorXr total_sum_epi = VectorXr::Zero(X.cols()); //del
	VectorXr total_sum_hipo = VectorXr::Zero(X.cols()); //del




	for(UInt triangle = 0; triangle < mesh_.num_elements(); ++triangle){

		MatrixXr X_cap = MatrixXr::Zero(Integrator::NNODES, 1);
		std::unique_ptr<Depth> depth_cap;

		Element<EL_NNODES, mydim, ndim> tri_activated = mesh_.getElement(triangle);

		MatrixXr sub_w = MatrixXr::Zero(EL_NNODES, X.cols());
		//bool contieneNan;

		for(Eigen::Index j=0; j < X.cols(); ++j){
			const VectorXr& x = X.col(j);
			Eigen::Matrix<Real, EL_NNODES, 1> sub_x;
			for(UInt i = 0; i < EL_NNODES; ++i){
				sub_x[i] = x[tri_activated[i].getId()];
			  	sub_w(i, j) = this->getWeights()(tri_activated[i].getId(), j);
			}

			Eigen::Matrix<Real, Integrator::NNODES, 1> x_cap =  MatrixXr::Zero(Integrator::NNODES, 1);
			for(UInt i=0; i<Integrator::NNODES; ++i){
				for(UInt j=0; j<EL_NNODES; ++j){
					if(isnan(sub_x[j]) && PsiQuad_(i,j) == 0)
						x_cap[i] += 0; 
					else
					x_cap[i] += sub_x[j]*PsiQuad_(i,j);
				}
			}
        
			X_cap.col(X_cap.cols()-1) = x_cap ; //sub_x.transpose() fill the matrix with the transformed functions (due to quadrature formula evaluation)

			if(x != X.col( X.cols()-1 ))
				X_cap.conservativeResize(X_cap.rows(), X_cap.cols()+1);
		}

		X_cap.conservativeResize(X_cap.rows(), X.cols());
		depth_cap = Depth_factory::createDepth(X_cap, d_tag);

		for(Eigen::Index j=0; j < X.cols(); ++j){
			Eigen::Matrix<Real, Integrator::NNODES, 1> weights = (PsiQuad_*sub_w.col(j)).array();

			//contieneNan = false;

			/* for(Eigen::Index i=0; i<EL_NNODES; ++i){
				if( isnan(X_cap(i,j)) ){ 
					contieneNan = true; 
				}
			} */

			Eigen::DiagonalMatrix<Real, Integrator::NNODES, Integrator::NNODES> weighdiag;
			weighdiag.diagonal() = weights;

			//if(true){


				/* VectorXr prova = VectorXr::Zero(Integrator::NNODES);
				prova = depth_cap->depth_epi(j);
				VectorXr prova2 = VectorXr::Zero(Integrator::NNODES);
				prova2 = depth_cap->depth_hipo(j); */

				if (d_tag == "MHRD"){
					total_sum_epi[j]  += (weighdiag*depth_cap->depth_epi(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
					total_sum_hipo[j]  += (weighdiag*depth_cap->depth_hipo(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
				}
				else{
					total_sum[j] += (weighdiag*depth_cap->depth(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0])) * tri_activated.getMeasure();
				}
			//}
			//else{
				//total_sum_epi[j]  += (weighdiag*depth_cap->depth_epi(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]));
				//total_sum_hipo[j] += (weighdiag*depth_cap->depth_hipo(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]));
			//}
		}
	}

		//Logica per considerare gli estremi
	   /* 	MatrixXr X_cap = MatrixXr::Zero(Integrator::NNODES, 1);
		std::unique_ptr<Depth> depth_cap;

		MatrixXr sub_w = MatrixXr::Zero(EL_NNODES, X.cols());
		bool contieneNan;
		Eigen::Vector2d id(0,X.rows()-1); 
					Rprintf("1"); 

		for(Eigen::Index j=0; j < X.cols(); ++j){
			const VectorXr& x = X.col(j);
			Rprintf("2");
			Eigen::Matrix<Real, EL_NNODES, 1> sub_x;
			for(UInt i = 0; i < EL_NNODES; ++i){
				Rprintf("3");
				sub_x[i] = x[id[i]];
				Rprintf("4");
			  	sub_w(i, j) = this->getWeights()(id[i], j);
			}
		Rprintf("5");
	
		Eigen::Matrix<Real, Integrator::NNODES, 1> x_cap = (PsiQuad_*sub_x).array();
		X_cap.col(X_cap.cols()-1) = sub_x.transpose() ;//x_cap; // fill the matrix with the transformed functions (due to quadrature formula evaluation)
			if(x != X.col( X.cols()-1 ))
				X_cap.conservativeResize(X_cap.rows(), X_cap.cols()+1);
		}
		X_cap.conservativeResize(X_cap.rows(), X.cols());


		Rprintf("dopo transformation");
		Rprintf("\n");
		for(auto i=0; i<X_cap.rows(); ++i)
			{
			for(auto jj=0; jj<X_cap.cols() ; ++jj)
				{
				Rprintf("%f ", X_cap(i,jj));
			}
		Rprintf("\n");
		}
		Rprintf("\n");
		depth_cap = Depth_factory::createDepth(X_cap, d_tag);
		for(Eigen::Index j=0; j < X.cols(); ++j){
		Eigen::Matrix<Real, Integrator::NNODES, 1> weights = (PsiQuad_*sub_w.col(j)).array();
		Eigen::DiagonalMatrix<Real, Integrator::NNODES, Integrator::NNODES> weighdiag;
		weighdiag.diagonal() = weights;

	
		total_sum_epi[j]  += (weighdiag*depth_cap->depth_epi(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))/X.rows() ;
		total_sum_hipo[j]  += (weighdiag*depth_cap->depth_hipo(j)).dot(EigenMap2WEIGHTS(&Integrator::WEIGHTS[0]))/X.rows();
		} */
		//fine logica

	if (d_tag == "MHRD"){
		for (Eigen::Index j=0; j < X.cols(); ++j)
		{
			total_sum[j] = std::min(total_sum_epi[j], total_sum_hipo[j]);
		}
	}

	return total_sum;
}

#endif /* __IFD_INTEGRATION_IMP_H__ */
