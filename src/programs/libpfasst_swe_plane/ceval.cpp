#include <iomanip>
#include <math.h>
#include <string>

#include <benchmarks_plane/SWEBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "SWE_Plane_TS_l_erk_n_erk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_ln_erk.hpp"
//#include "SWE_Plane_TS_l_rexi.hpp"

#include "ceval.hpp"
#include "cencap.hpp"

/**
 * Write file to data and return string of file name
 */
std::string write_file(
		       PlaneDataCtx  &i_ctx,
		       const PlaneData &i_planeData,
		       const char* i_name	///< name of output variable
		       )
{
  char buffer[1024];
  
  // get the pointer to the Simulation Variables object
  SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // create copy
  PlaneData planeData(i_planeData);

  // Write the data into the file
  const char* filename_template = simVars->misc.output_file_name.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);
  planeData.file_physical_saveData_ascii(buffer);
  
  return buffer;
}

/**
 *  Write the spectrum to file and return string of file name
 **/
std::string write_spectrum_to_file(
				   PlaneDataCtx &i_ctx,
				   const PlaneData &i_planeData,
				   const char* i_name
				   )
{
  std::cout << "Warning! write_spectrum_to_file on the plame has not beem implemented" << std::endl;
  // char buffer[1024];
  
  // // get the pointer to the Simulation Variables object
  // SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // // create copy
  // PlaneData planeData(i_planeData);

  // // Write the spectrum into the file
  // const char* filename_template = simVars->misc.output_file_name_prefix.c_str();
  // sprintf(buffer, 
  // 	  filename_template, 
  // 	  i_name, 
  // 	  simVars->timecontrol.current_timestep_size);
  // planeData.spectrum_file_write(buffer);
  
  // return buffer;
}

/** 
 *  Write the physical invariants to file
 **/
std::string write_physical_invariants_to_file(
					      PlaneDataCtx &i_ctx,
					      const char* i_name
					      )
{
  char buffer[1024];
  
  // get the pointer to the Simulation Variables object
  SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // get the vector of time and invariants
  const std::vector<double>& time               = i_ctx.get_time();
  const std::vector<double>& mass               = i_ctx.get_mass();
  const std::vector<double>& energy             = i_ctx.get_energy();
  const std::vector<double>& potentialEnstrophy = i_ctx.get_potential_enstrophy();

  // Write the spectrum into the file
  const char* filename_template = simVars->misc.output_file_name.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);

  std::ofstream file(buffer, std::ios_base::trunc);

  file << std::setprecision(20);

  for (unsigned int i = 0; i < time.size(); ++i) 
    file << time[i] << " " 
	 << mass[i] << " " 
	 << energy[i] << " " 
	 << potentialEnstrophy[i] 
	 << std::endl;

  file.close();
  
  return buffer;
}

/**
 * Write the residuals to file
 **/
std::string write_residuals_to_file(
				    PlaneDataCtx &i_ctx,
				    int i_rank,
				    const char* i_name
				    )
{
  char buffer[1024];
  
  // get the pointer to the Simulation Variables object
  SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // get the vector of residuals
  const std::vector<std::vector<double> >& residuals = i_ctx.get_residuals();

  // Write the spectrum into the file
  const char* filename_template = simVars->misc.output_file_name.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);

  std::ofstream file(buffer, std::ios_base::trunc);

  file << std::setprecision(20);

  for (unsigned int i = 0; i < residuals[i_rank].size(); ++i) 
    {
      file << i << " " 
	   << residuals[i_rank][i] << " " 
	   << std::endl;
      
    }
	  
  file.close();
  
  return buffer;
}



extern "C"
{
  // initialization of the variables (initial condition)
  void cinitial(
		PlaneDataCtx *i_ctx, 
		double i_t,
		double i_dt,
		PlaneDataVars *o_Y
		) 
  {
    PlaneData& h_Y    = o_Y->get_h();
    PlaneData& vort_Y = o_Y->get_vort();
    PlaneData& div_Y  = o_Y->get_div();

    // set the time stepping params
    i_ctx->setup_time_steps(
			    0,
			    i_dt
			    );

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());

    // get the configuration for this level
    PlaneDataConfig* data_config              = i_ctx->get_plane_data_config(o_Y->get_level());
    PlaneDataConfig* data_config_nodealiasing = i_ctx->get_plane_data_config_nodealiasing();

    // get the operator for this level
    PlaneOperators* op              = i_ctx->get_plane_operators(o_Y->get_level());
    PlaneOperators* op_nodealiasing = i_ctx->get_plane_operators_nodealiasing();

    // // // instantiate h, vort, and div without dealiasing to get the initial condition
    // PlaneData h_Y_nodealiasing(data_config_nodealiasing);
    // PlaneData vort_Y_nodealiasing(data_config_nodealiasing);
    // PlaneData div_Y_nodealiasing(data_config_nodealiasing);

    // h_Y_nodealiasing.physical_set_zero();
    // vort_Y_nodealiasing.physical_set_zero();
    // div_Y_nodealiasing.physical_set_zero();

    // // get the initial condition in h, vort, and div
    // PlaneBenchmarksCombined::setupInitialConditions(h_Y_nodealiasing, 
    // 						     vort_Y_nodealiasing, 
    // 						     div_Y_nodealiasing, 
    // 						     *simVars, 
    // 						     *op_nodealiasing);
    
    // h_Y.load_nodealiasing(h_Y_nodealiasing);
    // vort_Y.load_nodealiasing(vort_Y_nodealiasing);
    // div_Y.load_nodealiasing(div_Y_nodealiasing);


    SWEBenchmarksCombined::setupInitialConditions(h_Y, 
						  vort_Y, 
						  div_Y, 
						  *simVars, 
						  *op);

    // output the configuration
    simVars->outputConfig();
    
    double current_simulation_time = 0;
    int nsteps                     = 0;
   
    write_file(*i_ctx, h_Y,    "prog_h_init");
    write_file(*i_ctx, vort_Y, "prog_vort_init");
    write_file(*i_ctx, div_Y,  "prog_div_init");

    write_spectrum_to_file(*i_ctx, h_Y,    "init_spectrum_h");
    write_spectrum_to_file(*i_ctx, vort_Y, "init_spectrum_vort");
    write_spectrum_to_file(*i_ctx, div_Y,  "init_spectrum_div");

        
    // // get the timestepper 
    // SWE_Plane_TS_lg_irk_lc_n_erk* timestepper = i_ctx->get_lg_irk_lc_n_erk_timestepper();
    // //SWE_Plane_TS_ln_erk* timestepper = i_ctx->get_ln_erk_timestepper();
  
    // std::cout << "current_simulation_time = " << current_simulation_time 
    //  	      << " i_t = " << i_t 
    //  	      << std::endl;
    
    // std::cout << "i_dt = " << i_dt << std::endl;

    // // compute the reference solution (i.e., obtained with the reference time stepper)
    //  while (current_simulation_time < i_t)
    //    {
    //  	if (nsteps%120 == 0) 
    //  	  {
    //  	    std::cout << "current_simulation_time = "
    //  		      << current_simulation_time
    //  		      << std::endl;

    // // 	    // std::string name = "prog_h_ref_"+std::to_string(nsteps)+"_";
    // // 	    // write_file(*i_ctx, h_Y, name.c_str());
    // // 	    // name = "prog_vort_ref_"+std::to_string(nsteps)+"_";
    // // 	    // write_file(*i_ctx, vort_Y,name.c_str());
    // // 	    // name = "prog_div_ref_"+std::to_string(nsteps)+"_";
    // // 	    // write_file(*i_ctx, div_Y, name.c_str());

    // 	  }

    // 	// solve the implicit system using the Helmholtz solver
    // 	timestepper->run_timestep(
    // 				  h_Y,
    // 				  vort_Y,
    // 				  div_Y,
    // 				  i_dt,
    // 				  current_simulation_time
    // 				  );

    // 	if (simVars->sim.viscosity != 0)
    // 	  {
    // 	    double scalar = simVars->sim.viscosity*i_dt;
    // 	    double r      = simVars->sim.earth_radius;
	    
    // 	    h_Y  = h_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
    // 	    vort_Y = vort_Y.spectral_solve_helmholtz(1.0, -scalar, r);
    // 	    div_Y  = div_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
    // 	  }
	
    // 	current_simulation_time += i_dt;
    // 	nsteps += 1;
    //    }

    //  std::cout << "nsteps = " << nsteps << std::endl;

    //  write_spectrum_to_file(*i_ctx, h_Y,  "spectrum_h_ref");
    //  write_spectrum_to_file(*i_ctx, vort_Y, "spectrum_vort_ref");
    //  write_spectrum_to_file(*i_ctx, div_Y,  "spectrum_div_ref");

    //  write_file(*i_ctx, h_Y,  "prog_h_ref");
    //  write_file(*i_ctx, vort_Y, "prog_vort_ref");
    //  write_file(*i_ctx, div_Y,  "prog_div_ref");

    //  FatalError("stop the simulation");
    
  }

  // finalizes the time step when libpfasst is done 
  // currently does nothing else than outputting the solution
  void cfinal(
	      PlaneDataCtx *i_ctx, 
	      PlaneDataVars *i_Y,
	      int i_nnodes, 
	      int i_niters,
	      int i_rank,
	      int i_nprocs
	      ) 
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();
    
    const int& level_id = i_Y->get_level();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables()); 

    if (i_nprocs == 1)
      {
	std::string filename = "prog_h_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, h_Y, filename.c_str());
	
	filename = "prog_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, vort_Y, filename.c_str());
	
	filename = "prog_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, div_Y, filename.c_str());
	
	filename = "spectrum_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, vort_Y, filename.c_str());
	
	filename = "spectrum_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, div_Y, filename.c_str());

	filename = "spectrum_h_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, h_Y, filename.c_str());

	filename = "invariants_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_physical_invariants_to_file(*i_ctx, filename.c_str());

      }
    else if (i_rank == 0)
      {
	
	std::string filename = "prog_h_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, h_Y, filename.c_str());
	
	filename = "prog_vort_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, vort_Y, filename.c_str());
	
	filename = "prog_div_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_file(*i_ctx, div_Y, filename.c_str());
	
	filename = "spectrum_vort_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, vort_Y, filename.c_str());
	
	filename = "spectrum_div_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, div_Y, filename.c_str());

	filename = "spectrum_h_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_spectrum_to_file(*i_ctx, h_Y, filename.c_str());

	filename = "invariants_nprocs_"+std::to_string(i_nprocs)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_physical_invariants_to_file(*i_ctx, filename.c_str());
	
      }

    if (level_id == simVars->libpfasst.nlevels-1)
      {	
	std::string filename = "";
	if (i_nprocs == 1)
	  filename = "residuals_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	else
	  filename = "residuals_nprocs_"+std::to_string(i_nprocs)+"_current_proc_"+std::to_string(i_rank)+"_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
	write_residuals_to_file(*i_ctx, i_rank, filename.c_str());
      }
  }

  // evaluates the explicit (nonlinear) piece
  void ceval_f1(PlaneDataVars *i_Y,
		double i_t,
		PlaneDataCtx *i_ctx,
		PlaneDataVars *o_F1
		)
  {       
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    PlaneData& h_F1    = o_F1->get_h();
    PlaneData& vort_F1 = o_F1->get_vort();
    PlaneData& div_F1  = o_F1->get_div();

    // get the time step parameters
    SimulationVariables* simVars = i_ctx->get_simulation_variables();
    
    // return immediately if no nonlinear terms
    if (simVars->pde.use_linear_div == 1)
      {
	c_sweet_data_setval(o_F1, 0.0);
	return;
      }


    SWE_Plane_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());
    
    // compute the explicit nonlinear right-hand side
    timestepper->euler_timestep_update_nonlinear(
						 h_Y, 
						 vort_Y,
						 div_Y,
						 h_F1,
						 vort_F1,
						 div_F1, 
						 simVars->timecontrol.max_simulation_time
						 );

  }

  // evaluates the first implicit piece o_F2 = F2(i_Y)
  void ceval_f2 (
		 PlaneDataVars *i_Y, 
		 double i_t, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F2
		 ) 
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    PlaneData& h_F2    = o_F2->get_h();
    PlaneData& vort_F2 = o_F2->get_vort();
    PlaneData& div_F2  = o_F2->get_div();
    
    // initialize the right-hand side
    c_sweet_data_setval(o_F2, 0.0);

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();
  
    // get the explicit timestepper 
    SWE_Plane_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());
	
    // compute the linear right-hand side
    timestepper->euler_timestep_update_linear(
					      h_Y, 
					      vort_Y,
					      div_Y,
					      h_F2,
					      vort_F2,
					      div_F2,
					      simVars->timecontrol.max_simulation_time
					      );
  }

  // solves the first implicit system for io_Y
  // then updates o_F2 with the new value of F2(io_Y)
  void ccomp_f2 (
		 PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 PlaneDataVars *i_Rhs, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F2
		 ) 
  {
    PlaneData& h_Y    = io_Y->get_h();
    PlaneData& vort_Y = io_Y->get_vort();
    PlaneData& div_Y  = io_Y->get_div();

    const PlaneData& h_Rhs    = i_Rhs->get_h();
    const PlaneData& vort_Rhs = i_Rhs->get_vort();
    const PlaneData& div_Rhs  = i_Rhs->get_div();

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // first copy the rhs into the solution vector
    // this is needed to call the SWEET function run_ts_timestep
    h_Y    = h_Rhs;
    vort_Y = vort_Rhs;
    div_Y  = div_Rhs;
    

    if (simVars->libpfasst.use_rexi) 
      {
	// // get the rexi time stepper
	// SWE_Plane_TS_l_rexi* timestepper = i_ctx->get_l_rexi_timestepper(io_Y->get_level());

	// // solve the implicit system using the Helmholtz solver
	// timestepper->run_timestep(
	// 			  h_Y,
	// 			  vort_Y,
	// 			  div_Y,
	// 			  i_dt,
	// 			  simVars->timecontrol.max_simulation_time
	// 			  );
	
      }
    else 
      {
	// get the irk timestepper
	SWE_Plane_TS_l_irk* timestepper = i_ctx->get_l_irk_timestepper(io_Y->get_level());
	
	// solve the implicit system using the Helmholtz solver
	timestepper->run_timestep(
				  h_Y,
				  vort_Y,
				  div_Y,
				  i_dt,
				  simVars->timecontrol.max_simulation_time
				  );

      }

    // now recompute F2 with the new value of Y
    ceval_f2(
      	     io_Y, 
     	     i_t, 
     	     i_ctx, 
     	     o_F2
     	     );
    

    PlaneData& h_F2    = o_F2->get_h();
    PlaneData& vort_F2 = o_F2->get_vort();
    PlaneData& div_F2  = o_F2->get_div();

    phi_F2  = (phi_Y  - phi_Rhs)  / i_dt;
    vort_F2 = (vort_Y - vort_Rhs) / i_dt;
    div_F2  = (div_Y  - div_Rhs)  / i_dt;

  }

  // evaluates the second implicit piece o_F3 = F3(i_Y)
  // currently contains the implicit artificial viscosity term
  void ceval_f3 (
		 PlaneDataVars *i_Y, 
		 double i_t, 
		 int i_level,
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    const PlaneData& h_Y    = i_Y->get_h();
    const PlaneData& vort_Y = i_Y->get_vort();
    const PlaneData& div_Y  = i_Y->get_div();

    PlaneData& h_F3    = o_F3->get_h();
    PlaneData& vort_F3 = o_F3->get_vort();
    PlaneData& div_F3  = o_F3->get_div();

    // initialize F3 to zero in case no artificial viscosity
    c_sweet_data_setval(o_F3, 0.0);

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // no need to do anything if no artificial viscosity
    if (simVars->sim.viscosity == 0) 
      return;

    std::cout << "Warning: viscosity on the plane not implemented yet" << std::endl;

    // // get the parameters used to apply diffusion
    // const double r    = simVars->sim.earth_radius;
    // const double visc = simVars->sim.viscosity;

    // h_F3 = h_Y;
    // h_F3.spectral_update_lambda(
    // 				  [&](
    // 				      int n, int m,
    // 				      std::complex<double> &io_data
    // 				      )
    // 				  {
    // 				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
    // 				  }
    // 				  );

    // vort_F3 = vort_Y;
    // vort_F3.spectral_update_lambda(
    // 				  [&](
    // 				      int n, int m,
    // 				      std::complex<double> &io_data
    // 				      )
    // 				  {
    // 				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
    // 				  }
    // 				  );

    // div_F3 = div_Y;
    // div_F3.spectral_update_lambda(
    // 				  [&](
    // 				      int n, int m,
    // 				      std::complex<double> &io_data
    // 				      )
    // 				  {
    // 				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
    // 				  }
    // 				  );
  }

  // solves the second implicit system for io_Y
  // then updates o_F3 with the new value o_F3 = F3(io_Y)
  // currently solve the implicit system formed with artificial viscosity term
  void ccomp_f3 (
		 PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt,
		 int i_level,
		 PlaneDataVars *i_Rhs, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    PlaneData& h_Y    = io_Y->get_h();
    PlaneData& vort_Y = io_Y->get_vort();
    PlaneData& div_Y  = io_Y->get_div();

    PlaneData& h_Rhs    = i_Rhs->get_h();
    PlaneData& vort_Rhs = i_Rhs->get_vort();
    PlaneData& div_Rhs  = i_Rhs->get_div();

    // initialize F3 to zero in case no artificial viscosity
    c_sweet_data_setval(o_F3, 0.0);
            
    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // no need to do anything if no artificial viscosity
    if (simVars->sim.viscosity == 0)
      return;

    std::cout << "Warning: viscosity on the plame has not been implemented yet" << std::endl;

    // // get the parameters used to apply diffusion
    // const double scalar = simVars->sim.viscosity*i_dt;
    // const double r      = simVars->sim.earth_radius;

    // // solve (1-dt*visc*diff_op)*rhs = y
    // h_Y    = h_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r); 
    // vort_Y = vort_Rhs.spectral_solve_helmholtz(1.0, -scalar, r); 
    // div_Y  = div_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r); 
    
    // // now recompute F3 with the new value of Y
    
    // // ceval_f3(
    // //  	     io_Y, 
    // // 	     i_t, 
    // // 	     i_level,
    // // 	     i_ctx, 
    // // 	     o_F3
    // // 	     );
    
    // PlaneData& h_F3    = o_F3->get_h();
    // PlaneData& vort_F3 = o_F3->get_vort();
    // PlaneData& div_F3  = o_F3->get_div();

    // h_F3    = (h_Y  - h_Rhs)  / i_dt;
    // vort_F3 = (vort_Y - vort_Rhs) / i_dt;
    // div_F3  = (div_Y  - div_Rhs)  / i_dt;
  
  }

}
