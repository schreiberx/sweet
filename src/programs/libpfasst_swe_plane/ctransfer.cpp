#include "ctransfer.hpp"
#include "cencap.hpp"

extern "C"
{
  void c_sweet_data_restrict(
			     PlaneDataVars *io_y_coarse, 
			     PlaneDataVars *i_y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     double t) {
    c_sweet_data_copy(i_y_fine, 
		      io_y_coarse);
  }  

  void c_sweet_data_interpolate(
				PlaneDataVars *io_y_fine, 
				PlaneDataVars *i_y_coarse, 
				int i_level_fine,
				int i_level_fome,
				double t) {
    c_sweet_data_copy(i_y_coarse, 
		      io_y_fine);
  }
}

