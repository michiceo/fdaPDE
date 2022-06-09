#ifndef __IFD_FE_IMP_H__
#define __IFD_FE_IMP_H__

#include <map>

template<UInt ORDER, UInt mydim, UInt ndim>
void
FEIFD<ORDER, mydim, ndim>::calculateFunctionalBoxplot()
{

	VectorXr ifd_sorted = ifd;
	MatrixXr data_ = dataProblem_.data();
	int middle;
	size_t l = data_.cols();

	std::multimap<Real, VectorXr> ifdToFunction;
	std::map<int, VectorXr> outlier;

	for(Eigen::Index i = 0; i < l; ++i){

		ifdToFunction.insert( std::pair<Real, VectorXr>(ifd_sorted[i],data_.col(i)) ); // I have ordered pairs depth-function
	}

	// definisco la mediana come la funzione con ifd maggiore
	size_t c = ifdToFunction.count(ifdToFunction.crbegin()->first);

	median = VectorXr::Zero(data_.rows());

  // in case there are more function with the same max depth
    std::pair<std::multimap<Real, VectorXr>::iterator,std::multimap<Real, VectorXr>::iterator> mediane = ifdToFunction.equal_range(ifdToFunction.crbegin()->first);
	for(auto m = mediane.first; m != mediane.second; ++m){
		for(Eigen::Index i = 0; i < data_.rows(); ++i){
			median[i] += m->second[i]/c;
		}
	}

	// definisco i quartili come minimo e massimo della zona che contiene il 50% delle funzioni
	firstQuartile = median;
	thirdQuartile = median;
	// definisco l'iqr
	VectorXr interQuartileRange = VectorXr::Zero(data_.rows());
	// definisco i whisker
	VectorXr lowerWhiskerFittizio = VectorXr::Zero(data_.rows());
	VectorXr upperWhiskerFittizio = VectorXr::Zero(data_.rows());

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

	 		if(lowerWhisker[i] > it3->second[i] && outlier.find(middle) == outlier.end())
	 			lowerWhisker[i] = it3->second[i];
	 		if(upperWhisker[i] < it3->second[i] && outlier.find(middle) == outlier.end())
	 			upperWhisker[i] = it3->second[i];

				++middle;
	 	}
	}
}

#endif /* __IFD_FE_IMP_H__ */
