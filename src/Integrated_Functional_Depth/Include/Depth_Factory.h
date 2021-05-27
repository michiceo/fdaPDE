#ifndef __DEPTH_FACTORY_H__
#define __DEPTH_FACTORY_H__

#include <memory>
#include "Depth.h"

//! brief@ A Factory class: a class for the choice of the depth.

class Depth_factory {
public:
	static std::shared_ptr<Depth>
	createDepth(const MatrixXr& m, const std::string& d)
	{
		if(d == "MHRD")
			return std::make_shared<MHRD>(m);
		else
		{
			Rprintf("Unavailable depth choice - using Modified Half Region Depth");
			return std::make_shared<MHRD>(m);
		}
	}
};

#endif /* __DEPTH_FACTORY_H__ */

