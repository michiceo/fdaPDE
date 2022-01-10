#ifndef __IFD_DEPTH_FACTORY_H__
#define __IFD_DEPTH_FACTORY_H__

#include <memory>
#include "IFD_Depth.h"

//! brief@ A Factory class: a class for the choice of the depth.

class Depth_factory {
public:
	static std::shared_ptr<Depth>
	createDepth(const MatrixXr& m, const std::string& d)
	{
		if(d == "MHRD")
			return std::make_shared<MHRD>(m);
		else if(d == "MBD")
			return std::make_shared<MBD>(m);
		else
		{
			Rprintf("Unavailable depth choice - using Modified Band Depth");
			return std::make_shared<MBD>(m);
		}
	}
};

#endif /* __IFD_DEPTH_FACTORY_H__ */
