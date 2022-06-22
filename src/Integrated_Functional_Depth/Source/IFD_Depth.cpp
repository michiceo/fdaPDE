#include "../Include/IFD_Depth.h"
#include <math.h>

const UInt
Depth::isnan_vector(const VectorXr& v) const
{
  UInt number_nan=0;

  omp_set_num_threads(2);

  #pragma omp parallel for default(none) shared(v) reduction(+:number_nan)
  for(size_t i=0; i<v.size(); ++i)
    number_nan += isnan(v[i]);

  return number_nan;
}

const VectorXi
Depth::ranking(const VectorXr& v) const
{
  //omp_set_num_threads(2);

  VectorXr vsorted = v;
  VectorXi vrank   = VectorXi::Zero(v.size());

  std::sort(vsorted.data(), vsorted.data() + vsorted.size());

  std::map<double, int> ranks;
  int rank = 1;

  //#pragma omp parallel for default(none) shared(ranks, v, vsorted) reduction(+:rank)
  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    double element = vsorted[index];
    //#pragma omp atomic read
    ranks[element] = rank;
    rank++;
  }

  //#pragma omp parallel for default(none) shared(ranks, v, vrank)
  for(Eigen::Index index = 0; index < v.size(); index++)
  {
    double element = v[index];
    //#pragma omp atomic read
    vrank[index] = ranks[element];
  }
	return vrank;
}

MHRD::MHRD(const MatrixXr& m):
Depth(m)
{
}

const std::tuple<VectorXr, VectorXr, VectorXr>// VectorXr
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
	return std::make_tuple(mepi, mhipo, hrd);
}

const std::tuple<VectorXr, VectorXr, VectorXr>// VectorXr
MHRD::compute_depth(UInt k) const
{
	VectorXr mepi  = VectorXr::Zero(this->p_);
	VectorXr mhipo = VectorXr::Zero(this->p_);
	VectorXr hrd   = VectorXr::Zero(this->p_);

  omp_set_num_threads(2);

  // Real time1 = omp_get_wtime();

  # pragma omp parallel for default(none) shared(mepi, mhipo, hrd, p_, n_, m_, k)
  for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] > 0 ){
      #pragma omp atomic write
      mepi[j] = (double) (this->n_ - number_nan - rmat[k] + 1) / ((double) (this->n_ - number_nan)*this->p_);
      mhipo[j] = (double) (rmat[k]) / ((double) (this->n_ - number_nan)*this->p_);
		}
		else{
      #pragma omp atomic write
			mepi[j] = 0;
			mhipo[j] = 0;
		}
    #pragma omp atomic write
	  hrd[j] = std::min(mepi[j], mhipo[j]);
	}

  // Real time2 = omp_get_wtime();
  //
  // Real time12 = (time2-time1);
	// Rprintf("time12: %d\n", time12);

	return std::make_tuple(mepi, mhipo, hrd);
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

  omp_set_num_threads(2);

  // Real time1 = omp_get_wtime();

  //# pragma omp parallel for default(none) shared(mbd, p_, n_, m_, k)
	for(Eigen::Index j=0; j < this->p_; ++j){

    VectorXi rmat = ranking(m_.row(j));
    UInt number_nan = isnan_vector(m_.row(j));

    if(this->n_ - number_nan - rmat[k] >= 0 ){
      //#pragma omp atomic write
      mbd[j] = (double) ((this->n_ - number_nan - rmat[k])*(rmat[k] - 1) + (this->n_ - number_nan - 1)) / ((double) (this->n_ - number_nan)*(this->n_ - number_nan - 1)*this->p_/2);
		}
    else{
      //#pragma omp atomic write
			mbd[j] = 0;
		}
  }

  // Real time2 = omp_get_wtime();
  //
  // Real time12 = time2-time1;
	// Rprintf("time12: %lf\n", time12);

	return mbd;
}
