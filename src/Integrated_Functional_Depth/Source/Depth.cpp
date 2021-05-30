#include "../Include/Depth.h"

MHRD::MHRD(const MatrixXr& m):
Depth(m)
{
}

const VectorXi
MHRD::ranking(const VectorXr& v) const
{

	VectorXr vsorted = v;
  	VectorXi vrank   = VectorXi::Zero(v.size());

	std::sort(vsorted.data(), vsorted.data() + vsorted.size());

	for(Eigen::Index i=0; i<v.size(); ++i){
		auto it = std::find(v.data(), v.data() + v.size(), vsorted[i]);
		if(vrank[it - v.data()] != 0){
			it = std::find(it+1, v.data() + v.size(), vsorted[i]);
		}
		vrank[it - v.data()] = i+1;
	}
	return vrank;
}

const VectorXr
MHRD::compute_depth() const
{

	VectorXr mepi  = VectorXr::Zero(this->n_);
	VectorXr mhipo = VectorXr::Zero(this->n_);
	VectorXr hrd   = VectorXr::Zero(this->n_);

	for(Eigen::Index j=0; j < this->p_; ++j){

	  VectorXi rmat = ranking(m_.row(j));
	  
	  /*for(int k=0; k<rmat.size(); ++k){
	    Rprintf("%i \n", rmat[k]);
	  }*/
	    
		for(Eigen::Index i=0; i < this->n_; ++i){
			mepi[i]  = mepi[i]  + this->n_ - rmat[i];
			mhipo[i] = mhipo[i] + rmat[i]  - 1;
		}
	}

	std::for_each(mepi.data(), mepi.data() + mepi.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});
	std::for_each(mhipo.data(), mhipo.data() + mhipo.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});

	// Eigen::Index it = 0;
	/*std::for_each( hrd.data(), hrd.data() + hrd.size(), [&mepi, &mhipo, &it] (Real & r) {r = std::min(mepi[it++], mhipo[it++]);} );
  */
	
	for(Eigen::Index i=0; i < this->n_; ++i){
	  hrd[i] = std::min(mepi[i], mhipo[i]);
	}
	
	return hrd;

}

