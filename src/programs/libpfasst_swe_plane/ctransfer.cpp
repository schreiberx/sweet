#include "ctransfer.hpp"
#include "cencap.hpp"

#include <benchmarks_plane/SWEBenchmarksCombined.hpp>

extern "C"
{
  void c_sweet_data_restrict(
			     PlaneDataVars *io_Y_coarse, 
			     PlaneDataVars *i_Y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     PlaneDataCtx *i_ctx,
			     double i_t) 
  {
    const PlaneData& h_Y_fine    = i_Y_fine->get_h();
    const PlaneData& vort_Y_fine = i_Y_fine->get_vort();
    const PlaneData& div_Y_fine  = i_Y_fine->get_div();

    PlaneData& h_Y_coarse    = io_Y_coarse->get_h();
    PlaneData& vort_Y_coarse = io_Y_coarse->get_vort();
    PlaneData& div_Y_coarse  = io_Y_coarse->get_div();

    // restrict the fine variables to the coarse grid and copy into the coarse variables
    h_Y_coarse    = h_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
    vort_Y_coarse = vort_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
    div_Y_coarse  = div_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
 
  }  

  void c_sweet_data_interpolate(
				PlaneDataVars *io_Y_fine, 
				PlaneDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				PlaneDataCtx *i_ctx,
				double i_t) 
  {
    const PlaneData& h_Y_coarse    = i_Y_coarse->get_h();
    const PlaneData& vort_Y_coarse = i_Y_coarse->get_vort();
    const PlaneData& div_Y_coarse  = i_Y_coarse->get_div();
  
    PlaneData& h_Y_fine    = io_Y_fine->get_h();
    PlaneData& vort_Y_fine = io_Y_fine->get_vort();
    PlaneData& div_Y_fine  = io_Y_fine->get_div();

    // interpolate the coarse variables on the fine grid and copy into the coarse variables
    h_Y_fine    = h_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
    vort_Y_fine = vort_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
    div_Y_fine  = div_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));

  }
}

