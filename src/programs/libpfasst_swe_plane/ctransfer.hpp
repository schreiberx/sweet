#ifndef _CTRANSFER_HPP_
#define _CTRANSFER_HPP_

#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"

extern "C"
{
  void c_sweet_data_restrict(
			     PlaneDataVars *io_Y_coarse, 
			     PlaneDataVars *i_Y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     PlaneDataCtx *i_ctx,
			     double i_t
			     );

  void c_sweet_data_interpolate(
				PlaneDataVars *io_Y_fine, 
				PlaneDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				PlaneDataCtx *i_ctx,
				double i_t
				);

}

#endif
