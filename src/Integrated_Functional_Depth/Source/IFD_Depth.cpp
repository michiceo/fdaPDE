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

  std::map<Real, int> ranks;
  int rank = 1;

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    Real element = vsorted[index];
    ranks[element] = rank;
    rank++;
  }

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    Real element = v[index];
    vrank[index] = ranks[element];
  }
	return vrank;
}

MHRD::MHRD(const MatrixXr& m):
Depth(m)
{
}

const std::tuple<VectorXr, VectorXr, VectorXr>// VectorXr
MHRD::compute_depth(UInt k) const
{
  VectorXr mhrd = VectorXr::Zero(this->p_);
  VectorXr mepi = VectorXr::Zero(this->p_);
  VectorXr mhipo = VectorXr::Zero(this->p_);

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
	  mhrd[j] = std::min(mepi[j], mhipo[j]);
	}

	return std::make_tuple(mepi, mhipo, mhrd);
}



MBD::MBD(const MatrixXr& m):
Depth(m)
{
}

const VectorXr
MBD::compute_depth(UInt k) const
{
	VectorXr mbd = VectorXr::Zero(this->p_);

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
