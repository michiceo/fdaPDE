#ifndef __FE_INTEGRATED_FUNCTIONAL_DEPTH_IMP_H__
#define __FE_INTEGRATED_FUNCTIONAL_DEPTH_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
void
FEIFD<ORDER, mydim, ndim>::apply()
{
	ifd = dataProblem_.FEintegrate_depth(dataProblem_.data());
	depth = dataProblem_.getDepth();
}

#endif /* __FE_INTEGRATED_FUNCTIONAL_DEPTH_IMP_H__ */

