#ifndef _CTRANSFER_HPP_
#define _CTRANSFER_HPP_

#include "../libpfasst_interface/SphereDataVars.hpp"
#include "SphereDataCtxSDC.hpp"

extern "C"
{
  void c_sweet_data_restrict(
			     SphereDataVars *io_Y_coarse, 
			     SphereDataVars *i_Y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     SphereDataCtxSDC *i_ctx,
			     double i_t
			     );

  void c_sweet_data_interpolate(
				SphereDataVars *io_Y_fine, 
				SphereDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				SphereDataCtxSDC *i_ctx,
				double i_t
				);

}

#endif
