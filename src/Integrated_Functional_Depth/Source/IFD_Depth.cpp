#include "../Include/IFD_Depth.h"
#include <math.h>

const UInt
Depth::isnan_vector(const VectorXr& v) const
{
  UInt number_nan=0;

  for(size_t i=0; i<v.size(); ++i)
    number_nan += isnan(v[i]);

  return number_nan;
}


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
MHRD::compute_depth(UInt k) const
{
	VectorXr mepi  = VectorXr::Zero(this->p_);
	VectorXr mhipo = VectorXr::Zero(this->p_);
	VectorXr hrd   = VectorXr::Zero(this->p_);

	for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] > 0 ){
      mepi[j] = (double) (this->n_ - number_nan - rmat[k]) / ((double) (this->n_ - number_nan - 1));
      mhipo[j] = (double) (rmat[k] - 1) / ((double) (this->n_ - number_nan - 1));
		}
		else{
			mepi[j] = 0;
			mhipo[j] = 0;
		}
	  hrd[j] = std::min(mepi[j], mhipo[j]);
	}
	return hrd;
}

const VectorXr
MHRD::compute_depth() const
{
	VectorXr mepi  = VectorXr::Zero(this->n_);
	VectorXr mhipo = VectorXr::Zero(this->n_);
	VectorXr hrd   = VectorXr::Zero(this->n_);

	for(Eigen::Index j=0; j < this->p_; ++j){
		VectorXi rmat = ranking(m_.row(j));

	  UInt number_nan = isnan_vector(m_.row(j));

		for(Eigen::Index i=0; i < this->n_ - number_nan; ++i){
			mepi[i]  = mepi[i]  + (double)(this->n_ - number_nan - rmat[i])/((double) (this->n_ - number_nan - 1));
			mhipo[i] = mhipo[i] + (double)(rmat[i]  - 1)/((double) (this->n_ - number_nan - 1));
		}
	}

	// std::for_each(mepi.data(), mepi.data() + mepi.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});
	// std::for_each(mhipo.data(), mhipo.data() + mhipo.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});

	for(Eigen::Index i=0; i < this->n_; ++i){
		hrd[i] = std::min(mepi[i], mhipo[i]);
	}
	return hrd;
}
