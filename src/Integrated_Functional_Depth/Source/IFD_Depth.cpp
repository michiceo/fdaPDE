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

const VectorXi
Depth::ranking(const VectorXr& v) const
{
  VectorXr vsorted = v;
  VectorXi vrank   = VectorXi::Zero(v.size());

  std::sort(vsorted.data(), vsorted.data() + vsorted.size());

  std::map<double, int> ranks;
  int rank = 1;

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    double element = vsorted[index];
    ranks[element] = rank;
    rank++;
  }

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    double element = v[index];
    vrank[index] = ranks[element];
  }
	return vrank;
}

MHRD::MHRD(const MatrixXr& m):
Depth(m)
{
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
			mepi[i]  = mepi[i]  + (double)(this->n_ - number_nan - rmat[i] + 1)/((double) (this->n_ - number_nan)*this->p_);
			mhipo[i] = mhipo[i] + (double)(rmat[i])/((double) (this->n_ - number_nan)*this->p_);
		}
	}

	// std::for_each(mepi.data(), mepi.data() + mepi.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});
	// std::for_each(mhipo.data(), mhipo.data() + mhipo.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});

	for(Eigen::Index i=0; i < this->n_; ++i){
		hrd[i] = std::min(mepi[i], mhipo[i]);
	}
	return hrd;
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
      mepi[j] = (double) (this->n_ - number_nan - rmat[k] + 1) / ((double) (this->n_ - number_nan)*this->p_);
      mhipo[j] = (double) (rmat[k]) / ((double) (this->n_ - number_nan)*this->p_);
		}
		else{
			mepi[j] = 0;
			mhipo[j] = 0;
		}
	  hrd[j] = std::min(mepi[j], mhipo[j]);
	}
	return hrd;
}



MBD::MBD(const MatrixXr& m):
Depth(m)
{
}

const VectorXr
MBD::compute_depth() const
{
	VectorXr mbd  = VectorXr::Zero(this->n_);

	for(Eigen::Index j=0; j < this->p_; ++j){
		VectorXi rmat = ranking(m_.row(j));

	  UInt number_nan = isnan_vector(m_.row(j));

		for(Eigen::Index i=0; i < this->n_ - number_nan; ++i){
			mbd[i] += (double) ((this->n_ - number_nan - rmat[i])*(rmat[i] - 1) + (this->n_ - number_nan - 1)) / ((double) (this->n_ - number_nan)*(this->n_ - number_nan - 1)*this->p_/2);
    }
	}
	return mbd;
}

const VectorXr
MBD::compute_depth(UInt k) const
{
	VectorXr mbd   = VectorXr::Zero(this->p_);

	for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] >= 0 ){
      mbd[j] = (double) ((this->n_ - number_nan - rmat[k])*(rmat[k] - 1) + (this->n_ - number_nan - 1)) / ((double) (this->n_ - number_nan)*(this->n_ - number_nan - 1)*this->p_/2);
		}
    else{
			mbd[j] = 0;
		}
  }
	return mbd;
}
