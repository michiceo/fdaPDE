#include "../Include/IFD_Data.h"

IFDData::IFDData(SEXP Rdata, SEXP Rorder, SEXP Rweights)
{
	setData(Rdata);
	order_ = INTEGER(Rorder)[0];
	setWeights(Rweights);
}

IFDData::IFDData(const MatrixXr & data, const UInt & order, const VectorXr & weights):
		data_(data), order_(order), weights_(weights)
{
}

void
IFDData::setData(SEXP Rdata)
{
  const RNumericMatrix data(Rdata);
  UInt n_=data.nrows(), p_=data.ncols();
  if(n_>0){
	  data_.resize(n_, p_);
	  for(auto i=0; i<n_; ++i)
	  {
	  	for(auto j=0; j<p_ ; ++j)
	  	{
	  		data_(i,j)=REAL(Rdata)[i+n_*j];
	  	}
	  }
  }
}

void
IFDData::setWeights(SEXP Rweights)
{
  UInt dimc = Rf_length(Rweights);
  weights_.resize(dimc);
  for(UInt i=0; i<dimc; i++)
  {
      weights_[i] = REAL(Rweights)[i];
  }
}

void
IFDData::printData(std::ostream & out) const
{
	for(auto i=0;i<data_.rows(); i++)
	{
		for(auto j=0;j<data_.cols();j++)
		{
		out<<data_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

