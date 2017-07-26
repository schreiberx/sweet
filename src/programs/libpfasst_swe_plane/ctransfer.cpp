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
    PlaneData& h_Y_fine = i_Y_fine->get_h();
    const PlaneData& u_Y_fine = i_Y_fine->get_u();
    const PlaneData& v_Y_fine = i_Y_fine->get_v();

    PlaneData& h_Y_coarse = io_Y_coarse->get_h();
    PlaneData& u_Y_coarse = io_Y_coarse->get_u();
    PlaneData& v_Y_coarse = io_Y_coarse->get_v();

    // //////////////////////////////////////////////

    // // Try this with and without dealiasing

    // SimulationVariables* simVars(i_ctx->get_simulation_variables());
    
    // h_Y_coarse.physical_update_lambda_array_indices(
    // 						 [&](int i, int j, double &o_data)
    // 						 {
    // 						   o_data = (i+1.0)+(j+3.0)*j;
						   
    // 						   /*
    // 						   double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
    // 						   o_data *= std::cos(((double)x)*(double)fx*M_PI*2.0/(double)simVars->disc.res_spectral[1]);
    // 						   */
    // 						 }
    // 						 );
    // // print the original data
    // h_Y_coarse.print_physicalArrayData();

    // // create a new PlaneData object and print its data
    // PlaneData tmp(i_ctx->get_plane_data_config(i_level_coarse));
    // tmp = h_Y_coarse;
    // tmp.print_physicalArrayData();
		
    // // convert to spectral data and then back to spectral
    // tmp.request_data_spectral();
    // tmp.request_data_physical();

    // // manually compute and output the difference between the two PlaneData objects
    // for (std::size_t i = 0; i < i_ctx->get_plane_data_config(i_level_coarse)->physical_array_data_number_of_elements; i++)
    //   std::cout << tmp.physical_space_data[i] - h_Y_coarse.physical_space_data[i] << std::endl;    
    
    // // interpolate and then restrict 
    // h_Y_fine = h_Y_coarse.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_fine));
    // h_Y_coarse = h_Y_fine.spectral_returnWithDifferentModes(i_ctx->get_plane_data_config(i_level_coarse));

    // // print the result
    // h_Y_coarse.print_physicalArrayData();
    
    // // compute the error
    // tmp = h_Y_coarse-tmp;

    // // compute the norm
    // const double error = tmp.reduce_maxAbs();
    
    // // print the norm and return
    // std::cout << "error = " << error << std::endl;
    // exit(0);
    
    // //////////////////////////////////////////////

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

