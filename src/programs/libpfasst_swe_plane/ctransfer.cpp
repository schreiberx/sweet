#include <sweet/plane/PlaneData.hpp>
#include "PlaneDataCtx.hpp"

extern "C"
{
  void c_sweet_data_restrict(
			     PlaneData *Y_G, 
			     PlaneData *Y_F, 
			     double t) {
    // not implemented
  }  

  void c_sweet_data_interpolate(
				PlaneData *Y_F, 
				PlaneData *Y_G, 
				double t) {
    // not implemented
  }
}

