#include "../Include/IFD_Depth.h"
#include <math.h>
#include <R.h>


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

    auto cmp=[](const double &d1, const double &d2){
    if (isnan(d1)) return false;
    if (isnan(d2)) return true;
    return d1<d2;
    };  

  std::sort(vsorted.data(), vsorted.data() + vsorted.size(), cmp);
  std::map<Real, int> ranks;
  int rank = 0;

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    Real element = vsorted[index];
    if (!isnan(element)){
      ranks[element] = rank;
      rank++;
    }
    
  }

  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    Real element = v[index];
    if (!isnan(element)){
      vrank[index] = ranks[element];
    }
    else{
      vrank[index] = v.size();
    }
    //Rprintf("- %f %i %f - ", v[index] , ranks[element], v[index]);
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

    if(this->n_ - number_nan - rmat[k] > 0 ){  // ci va = ?
      mepi[j] = (double) (this->n_ - number_nan - rmat[k] - 1 )/(double) (this->n_ - number_nan) ;// / /*(*/(double) (this->n_ - number_nan)/* *this->p_) */ ;
      mhipo[j] = (double) (rmat[k] )/(double) (this->n_ - number_nan) ;// / /*(*/(double) (this->n_ - number_nan)/* *this->p_) */ ;
		}
		else{//qua se sono NA
			mepi[j] = 0;
			mhipo[j] = 0;
		}
	  mhrd[j] = std::min(mepi[j], mhipo[j]);

  }
	return std::make_tuple(mepi, mhipo, mhrd);
}


MBD::MBD(const MatrixXr& m):
Depth(m){}

const VectorXr
MBD::compute_depth(UInt k) const
{
	VectorXr mbd = VectorXr::Zero(this->p_);

	for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] > 0 ){
      mbd[j] = (double) 2*(this->n_ - number_nan - rmat[k] - 1)*(rmat[k] + 1)/(double) ((this->n_ - number_nan)*(this->n_ - number_nan));
      //mbd[j] = (double) ((this->n_ - number_nan - rmat[k])*(rmat[k] - 1) + (this->n_ - number_nan - 1)) / ((double) (this->n_ - number_nan)*(this->n_ - number_nan - 1)*this->p_/2);
    }
    else{
			mbd[j] = 0;
		}
  }
	return mbd;
}

FMD::FMD(const MatrixXr& m):
Depth(m){}

const VectorXr
FMD::compute_depth(UInt k) const
{
	VectorXr fmd = VectorXr::Zero(this->p_);

	for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] > 0 ){
      fmd[j] = 1 - abs((double) (1-2*(rmat[k] + 1))/(double) (2*(this->n_ - number_nan)));
      //mbd[j] = (double) ((this->n_ - number_nan - rmat[k])*(rmat[k] - 1) + (this->n_ - number_nan - 1)) / ((double) (this->n_ - number_nan)*(this->n_ - number_nan - 1)*this->p_/2);
    }
    else{
			fmd[j] = 0;
		}
  }
	return fmd;
}
