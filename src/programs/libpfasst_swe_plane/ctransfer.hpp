#ifndef _CTRANSFER_HPP_
#define _CTRANSFER_HPP_

#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"

extern "C"
{
  void c_sweet_data_restrict(
			     PlaneDataVars *io_y_coarse, 
			     PlaneDataVars *i_y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     double t
			     );

  void c_sweet_data_interpolate(
				PlaneDataVars *io_y_fine, 
				PlaneDataVars *i_y_coarse, 
				int i_level_fine,
				int i_level_fome,
				double t
				);

}

#endif
