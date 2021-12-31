#include "ctransfer.hpp"

#include "../swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"

#include "cencap.hpp"


extern "C"
{
  void c_sweet_data_restrict(
			     SphereDataVars *io_Y_coarse, 
			     SphereDataVars *i_Y_fine, 
			     int i_level_coarse,
			     int i_level_fine, 
			     SphereDataCtxSDC *i_ctx,
			     double i_t) 
  {
    const SphereData_Spectral& phi_pert_Y_fine  = i_Y_fine->get_phi_pert();
    const SphereData_Spectral& vrt_Y_fine = i_Y_fine->get_vrt();
    const SphereData_Spectral& div_Y_fine  = i_Y_fine->get_div();

    SphereData_Spectral& phi_pert_Y_coarse  = io_Y_coarse->get_phi_pert();
    SphereData_Spectral& vrt_Y_coarse = io_Y_coarse->get_vrt();
    SphereData_Spectral& div_Y_coarse  = io_Y_coarse->get_div();

    // restrict the fine variables to the coarse grid and copy into the coarse variables
    phi_pert_Y_coarse  = phi_pert_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());
    vrt_Y_coarse = vrt_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());
    div_Y_coarse  = div_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());
 
  }  

  void c_sweet_data_interpolate(
				SphereDataVars *io_Y_fine, 
				SphereDataVars *i_Y_coarse, 
				int i_level_fine,
				int i_level_coarse,
				SphereDataCtxSDC *i_ctx,
				double i_t) 
  {
    const SphereData_Spectral& phi_pert_Y_coarse  = i_Y_coarse->get_phi_pert();
    const SphereData_Spectral& vrt_Y_coarse = i_Y_coarse->get_vrt();
    const SphereData_Spectral& div_Y_coarse  = i_Y_coarse->get_div();
  
    SphereData_Spectral& phi_pert_Y_fine  = io_Y_fine->get_phi_pert();
    SphereData_Spectral& vrt_Y_fine = io_Y_fine->get_vrt();
    SphereData_Spectral& div_Y_fine  = io_Y_fine->get_div();

    // interpolate the coarse variables on the fine grid and copy into the coarse variables
    phi_pert_Y_fine  = phi_pert_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());
    vrt_Y_fine = vrt_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());
    div_Y_fine  = div_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_sphere_data_config());

  }
}

