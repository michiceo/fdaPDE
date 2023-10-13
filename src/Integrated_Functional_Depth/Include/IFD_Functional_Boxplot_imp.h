#ifndef __IFD_FUNCTIONAL_BOXPLOT_IMP_H__
#define __IFD_FUNCTIONAL_BOXPLOT_IMP_H__

#include <map>

typedef std::multimap<Real, VectorXr>::iterator multimap_iterator;

template<UInt ORDER, UInt mydim, UInt ndim>
void
FunctionalBoxplot<ORDER, mydim, ndim>::calculateFunctionalBoxplot()
{

	const VectorXr& ifd_sorted = ifd;
	const MatrixXr& data_ = depthIntegration_.data();
	int middle;
	int l = data_.cols();

	std::multimap<Real, VectorXr> ifdToFunction;
	std::map<int, VectorXr> outlier;

	for(Eigen::Index i = 0; i < l; ++i){

		ifdToFunction.insert( std::pair<Real, VectorXr>(ifd_sorted[i],data_.col(i)) ); // we ordered pairs depth-function
	}

	// median is defined as the function with max depth value
	size_t c = ifdToFunction.count(ifdToFunction.crbegin()->first);

	median = VectorXr::Zero(data_.rows());

  // in case there are more function with the same max depth:
  std::pair<multimap_iterator, multimap_iterator> mediane = ifdToFunction.equal_range(ifdToFunction.crbegin()->first);

	for(auto m = mediane.first; m != mediane.second; ++m){
		for(Eigen::Index i = 0; i < data_.rows(); ++i){
			median[i] += m->second[i]/c;
		}
	}

	// quartiles are min and max of the area containing 50% of data
	firstQuartile = median;
	thirdQuartile = median;
	// IQR
	VectorXr interQuartileRange = VectorXr::Zero(data_.rows());
	// whisker
	VectorXr lowerWhiskerFittizio = VectorXr::Zero(data_.rows());
	VectorXr upperWhiskerFittizio = VectorXr::Zero(data_.rows());

	//middle = 0;

	// std::for_each(mepi.data(), mepi.data() + mepi.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});
	// std::for_each(mhipo.data(), mhipo.data() + mhipo.size(), [this] (Real & r) {r /= (double) (this->n_ - 1);});

	for(Eigen::Index i = 0; i < data_.rows(); ++i){

		middle = 0;

		for (auto it1 = ifdToFunction.rbegin(); middle < l/2 && it1 != ifdToFunction.rend(); ++it1){

	 		if(firstQuartile[i] > it1->second[i])
	 			firstQuartile[i] = it1->second[i];
	 		if(thirdQuartile[i] < it1->second[i])
	 			thirdQuartile[i] = it1->second[i];
	 		interQuartileRange[i] = thirdQuartile[i] - firstQuartile[i];
	 		upperWhiskerFittizio[i] = thirdQuartile[i] + 1.5*interQuartileRange[i];
	 		lowerWhiskerFittizio[i] = firstQuartile[i] - 1.5*interQuartileRange[i];

			++middle;

	 	}
	}

	for(Eigen::Index i = 0; i < data_.rows(); ++i){
		middle = 0;
	 	for (auto it2 = ifdToFunction.begin(); middle < l/2 && it2 != ifdToFunction.end(); ++it2){

	 		if(lowerWhiskerFittizio[i] > it2->second[i]){
	 			outlier.insert( std::pair<int, VectorXr>(middle,it2->second));
	 		}
	 		if(upperWhiskerFittizio[i] < it2->second[i]){
	 			outlier.insert( std::pair<int, VectorXr>(middle,it2->second));
	 		}
			++middle;
	 	}
	}

	lowerWhisker = firstQuartile;
	upperWhisker = thirdQuartile;

	for(Eigen::Index i = 0; i < data_.rows(); ++i){
		middle = 0;
	 	for (auto it3 = ifdToFunction.begin(); middle < l/2 && it3 != ifdToFunction.end(); ++it3){

	 		if(lowerWhisker[i] > it3->second[i] && outlier.find(middle) == outlier.end()){
	 			lowerWhisker[i] = it3->second[i];
			}
	 		if(upperWhisker[i] < it3->second[i] && outlier.find(middle) == outlier.end()){
	 			upperWhisker[i] = it3->second[i];
			}
				++middle;
	 	}
	}

	//Per capire quale funzione è maggiore

	double differece;
	signDepth = VectorXr::Zero(data_.cols());


	for(Eigen::Index i = 0; i < l; ++i){
		differece = 0;

		for(Eigen::Index j = 0; j < data_.rows(); ++j){
			differece = differece + data_(j,i) - median[j];
	 	}

		if (differece >= 0){
			signDepth[i] = 1;		
		}
		else{
			signDepth[i] = 0;		
		}
		
	}

}

	#endif /* __IFD_FUNCTIONAL_BOXPLOT_IMP_H__ */
