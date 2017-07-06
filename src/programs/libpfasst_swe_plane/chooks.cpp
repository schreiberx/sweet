#include <sweet/plane/PlaneData.hpp>
#include "PlaneDataCtx.hpp"

extern "C"
{
  void cecho_error(PlaneData* sd, 
	    int step) {
    // not implemented
  }

  void cecho_residual(PlaneData* sd, 
                      int step, 
	              int iter) {
    // not impllmented
  }

}
