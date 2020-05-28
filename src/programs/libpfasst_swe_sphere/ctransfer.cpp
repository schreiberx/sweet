#include "ctransfer.hpp"

#include <benchmarks_sphere/SWESphereBenchmarks.hpp>

#include "cencap.hpp"


extern "C"
{
  void c_sweet_data_restrict(
			     SphereDataVars *io_Y_coarse, 
			     SphereDataVars *i_Y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     SphereDataCtx *i_ctx,
			     double i_t) 
  {
    const SphereData_Spectral& phi_Y_fine  = i_Y_fine->get_phi();
    const SphereData_Spectral& vort_Y_fine = i_Y_fine->get_vort();
    const SphereData_Spectral& div_Y_fine  = i_Y_fine->get_div();

    SphereData_Spectral& phi_Y_coarse  = io_Y_coarse->get_phi();
    SphereData_Spectral& vort_Y_coarse = io_Y_coarse->get_vort();
    SphereData_Spectral& div_Y_coarse  = io_Y_coarse->get_div();

    // restrict the fine variables to the coarse grid and copy into the coarse variables
    phi_Y_coarse  = phi_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_coarse));
    vort_Y_coarse = vort_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_coarse));
    div_Y_coarse  = div_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_coarse));
 
  }  

  void c_sweet_data_interpolate(
				SphereDataVars *io_Y_fine, 
				SphereDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				SphereDataCtx *i_ctx,
				double i_t) 
  {
    const SphereData_Spectral& phi_Y_coarse  = i_Y_coarse->get_phi();
    const SphereData_Spectral& vort_Y_coarse = i_Y_coarse->get_vort();
    const SphereData_Spectral& div_Y_coarse  = i_Y_coarse->get_div();
  
    SphereData_Spectral& phi_Y_fine  = io_Y_fine->get_phi();
    SphereData_Spectral& vort_Y_fine = io_Y_fine->get_vort();
    SphereData_Spectral& div_Y_fine  = io_Y_fine->get_div();

    // interpolate the coarse variables on the fine grid and copy into the coarse variables
    phi_Y_fine  = phi_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_fine));
    vort_Y_fine = vort_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_fine));
    div_Y_fine  = div_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config(i_level_fine));

  }
}

