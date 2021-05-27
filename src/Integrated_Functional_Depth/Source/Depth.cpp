#include "../Include/Depth.h"

MHRD::MHRD(const MatrixXr& m):
Depth(m)
{
}

const std::vector<Eigen::Index>
MHRD::ranking(const VectorXr& v){

	VectorXr vsorted = v;
	VectorXi vrank   = VectorXi::Zero(v.size());

	std::sort(vsorted.data(), vsorted.data() + vsorted.size());

	for(Eigen::Index i=0; i<v.size(); ++i){
		auto it = std::find(v.data(), v.data() + v.size(), vsorted[i]);
		if(vrank[it] != 0){
			it = std::find(it, v.data() + v.size(), vsorted[i]);
		}
		vrank[it] = i+1;
	}
	return vrank;
}

const VectorXr
MHRD::compute_depth() const
{

	VectorXr mepi  = VectorXr::Zero(n_);
	VectorXr mhipo = VectorXr::Zero(n_);
	VectorXr hrd   = VectorXr::Zero(n_);

	for (auto row : dataProblem_.data().rowwise()){

		std::vector<Eigen::Index> rmat = ranking(row);

		for(Eigen::Index i=0; i<n_; ++i){
			mepi[i]  = mepi[i]  + n_      - rmat[i];
			mhipo[i] = mhipo[i] + rmat[i] - 1;
		}
	}

	std::for_each(mepi.data(), mepi.data() + mepi.size(), [&n_] (Real & r) {r /= (double) n_;});
	std::for_each(mhipo.data(), mhipo.data() + mhipo.size(), [&n_] (Real & r) {r /= (double) n_;});

	Eigen::Index it = 0;
	std::for_each( hrd.data(), hrd.data() + hrd.size(), [&mepi, &mhipo, &it] (Real & r) {r = std::min(mepi[it++], mhipo[it++]);} );

	return hrd;

}

