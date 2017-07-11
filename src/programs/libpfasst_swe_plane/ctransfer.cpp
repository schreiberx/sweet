#include "ctransfer.hpp"
#include "cencap.hpp"

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
    const PlaneData& h_Y_fine = i_Y_fine->get_h();
    const PlaneData& u_Y_fine = i_Y_fine->get_u();
    const PlaneData& v_Y_fine = i_Y_fine->get_v();

    PlaneData& h_Y_coarse = io_Y_coarse->get_h();
    PlaneData& u_Y_coarse = io_Y_coarse->get_u();
    PlaneData& v_Y_coarse = io_Y_coarse->get_v();

    // restrict the fine variables to the coarse grid and copy into the coarse variables
    h_Y_coarse = h_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
    u_Y_coarse = u_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
    v_Y_coarse = v_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));
  }  

  void c_sweet_data_interpolate(
				PlaneDataVars *io_Y_fine, 
				PlaneDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				PlaneDataCtx *i_ctx,
				double i_t) 
  {
    const PlaneData& h_Y_coarse = i_Y_coarse->get_h();
    const PlaneData& u_Y_coarse = i_Y_coarse->get_u();
    const PlaneData& v_Y_coarse = i_Y_coarse->get_v();
  
    PlaneData& h_Y_fine = io_Y_fine->get_h();
    PlaneData& u_Y_fine = io_Y_fine->get_u();
    PlaneData& v_Y_fine = io_Y_fine->get_v();

    // interpolate the coarse variables on the fine grid and copy into the coarse variables
    h_Y_fine = h_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
    u_Y_fine = u_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
    v_Y_fine = v_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
  }
}

