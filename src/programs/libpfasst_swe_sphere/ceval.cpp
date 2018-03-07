#include <iomanip>
#include <math.h>
#include <string>

#include <benchmarks_sphere/SphereBenchmarksCombined.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>

#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_t_erk.hpp"
#include "SWE_Sphere_TS_lg_irk_lc_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_irk_n_erk_ver01.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_ln_erk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"

#include "ceval.hpp"
#include "cencap.hpp"

/**
 * Write file to data and return string of file name
 */
std::string write_file(
		       SphereDataCtx  &i_ctx,
		       const SphereData &i_sphereData,
		       const char* i_name	///< name of output variable
		       )
{
  char buffer[1024];
  
  // get the pointer to the Simulation Variables object
  SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // Write the data into the file
  const char* filename_template = simVars->misc.output_file_name_prefix.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);
  i_sphereData.physical_file_write(buffer);
  
  return buffer;
}

/**
 *  Write the spectrum to file and return string of file name
 **/
std::string write_spectrum_to_file(
				   SphereDataCtx &i_ctx,
				   const SphereData &i_sphereData,
				   const char* i_name
				   )
{
  char buffer[1024];
  
  // get the pointer to the Simulation Variables object
  SimulationVariables* simVars = i_ctx.get_simulation_variables();

  // Write the spectrum into the file
  const char* filename_template = simVars->misc.output_file_name_prefix.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);
  i_sphereData.spectrum_file_write(buffer);
  
  return buffer;
}

/** 
 *  Write the physical invariants to file
 **/
std::string write_physical_invariants_to_file(
					      SphereDataCtx &i_ctx,
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
  const char* filename_template = simVars->misc.output_file_name_prefix.c_str();
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


extern "C"
{
  // initialization of the variables (initial condition)
  void cinitial(
		SphereDataCtx *i_ctx, 
		double i_t,
		double i_dt,
		SphereDataVars *o_Y
		) 
  {
    SphereData& phi_Y  = o_Y->get_phi();
    SphereData& vort_Y = o_Y->get_vort();
    SphereData& div_Y  = o_Y->get_div();

    // set the time stepping params
    i_ctx->setup_time_steps(
			    i_t,
			    i_dt
			    );

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());

    if (simVars->sim.use_topography)
      write_file(*i_ctx, simVars->sim.h_topo,  "prog_h_topo");

    // initialize the variables
    phi_Y.physical_set_zero();
    vort_Y.physical_set_zero();
    div_Y.physical_set_zero();

    // get the configuration for this level
    SphereDataConfig* data_config              = i_ctx->get_sphere_data_config(o_Y->get_level());
    SphereDataConfig* data_config_nodealiasing = i_ctx->get_sphere_data_config_nodealiasing();

    // get the operator for this level
    SphereOperators* op              = i_ctx->get_sphere_operators(o_Y->get_level());
    SphereOperators* op_nodealiasing = i_ctx->get_sphere_operators_nodealiasing();

    // instantiate h, u, and v to get the initial condition
    SphereData phi_Y_nodealiasing(data_config_nodealiasing);
    SphereData vort_Y_nodealiasing(data_config_nodealiasing);
    SphereData div_Y_nodealiasing(data_config_nodealiasing);

    // get the initial condition in h, u, and v
    SphereBenchmarksCombined::setupInitialConditions(phi_Y_nodealiasing, 
						     vort_Y_nodealiasing, 
						     div_Y_nodealiasing, 
						     *simVars, 
						     *op_nodealiasing);
    
    phi_Y.load_nodealiasing(phi_Y_nodealiasing);
    vort_Y.load_nodealiasing(vort_Y_nodealiasing);
    div_Y.load_nodealiasing(div_Y_nodealiasing);
    
    // output the configuration
    simVars->outputConfig();
    
    double current_simulation_time = 0;
    int nsteps                     = 0;
   
    write_file(*i_ctx, phi_Y,  "prog_phi_init");
    write_file(*i_ctx, vort_Y, "prog_vort_init");
    write_file(*i_ctx, div_Y,  "prog_div_init");

    write_spectrum_to_file(*i_ctx, phi_Y,  "init_spectrum_phi");
    write_spectrum_to_file(*i_ctx, vort_Y, "init_spectrum_vort");
    write_spectrum_to_file(*i_ctx, div_Y,  "init_spectrum_div");

        
    // get the timestepper 
    //SWE_Sphere_TS_lg_irk_lc_n_erk* timestepper = i_ctx->get_lg_irk_lc_n_erk_timestepper();
    SWE_Sphere_TS_ln_erk* timestepper = i_ctx->get_ln_erk_timestepper();
  
    std::cout << "current_simulation_time = " << current_simulation_time 
    	      << " i_t = " << i_t 
    	      << std::endl;
    
    std::cout << "i_dt = " << i_dt << std::endl;

    // compute the reference solution (i.e., obtained with the reference time stepper)
    while (current_simulation_time < i_t)
      {
    	if (nsteps%100 == 0)
    	  std::cout << "current_simulation_time = "
    		    << current_simulation_time
    		    << std::endl;

    	// solve the implicit system using the Helmholtz solver
    	timestepper->run_timestep(
    				  phi_Y,
    				  vort_Y,
    				  div_Y,
    				  i_dt,
    				  simVars->timecontrol.max_simulation_time
    				  );

    	if (simVars->sim.viscosity != 0)
    	  {
    	    double scalar = simVars->sim.viscosity*i_dt;
    	    double r      = simVars->sim.earth_radius;
	    
    	    phi_Y  = phi_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
    	    vort_Y = vort_Y.spectral_solve_helmholtz(1.0, -scalar, r);
    	    div_Y  = div_Y.spectral_solve_helmholtz(1.0,  -scalar, r);
    	  }
	
    	current_simulation_time += i_dt;
    	nsteps += 1;
      }

    std::cout << "nsteps = " << nsteps << std::endl;

    write_spectrum_to_file(*i_ctx, phi_Y,  "spectrum_phi_ref");
    write_spectrum_to_file(*i_ctx, vort_Y, "spectrum_vort_ref");
    write_spectrum_to_file(*i_ctx, div_Y,  "spectrum_div_ref");


    write_file(*i_ctx, phi_Y,  "prog_phi_ref");
    write_file(*i_ctx, vort_Y, "prog_vort_ref");
    write_file(*i_ctx, div_Y,  "prog_div_ref");

    FatalError("stop the simulation");
    
  }

  // finalizes the time step when libpfasst is done 
  // currently does nothing else than outputting the solution
  void cfinal(
	      SphereDataCtx *i_ctx, 
	      SphereDataVars *i_Y,
	      int i_nnodes, 
	      int i_niters
	      ) 
  {
    const SphereData& phi_Y  = i_Y->get_phi();
    const SphereData& vort_Y = i_Y->get_vort();
    const SphereData& div_Y  = i_Y->get_div();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables()); 
    
    std::string filename = "prog_phi_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, phi_Y, filename.c_str());
    
    filename = "prog_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, vort_Y, filename.c_str());

    filename = "prog_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, div_Y, filename.c_str());
    
    filename = "spectrum_vort_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_spectrum_to_file(*i_ctx, vort_Y, filename.c_str());

    filename = "spectrum_div_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_spectrum_to_file(*i_ctx, div_Y, filename.c_str());

    filename = "spectrum_phi_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_spectrum_to_file(*i_ctx, phi_Y, filename.c_str());

    filename = "invariants_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_physical_invariants_to_file(*i_ctx, filename.c_str());
  
  }

  // evaluates the explicit (nonlinear) piece
  void ceval_f1(SphereDataVars *i_Y,
		double i_t,
		SphereDataCtx *i_ctx,
		SphereDataVars *o_F1
		)
  {       
    const SphereData& phi_Y  = i_Y->get_phi();
    const SphereData& vort_Y = i_Y->get_vort();
    const SphereData& div_Y  = i_Y->get_div();

    SphereData& phi_F1  = o_F1->get_phi();
    SphereData& vort_F1 = o_F1->get_vort();
    SphereData& div_F1  = o_F1->get_div();

    // get the time step parameters
    SimulationVariables* simVars = i_ctx->get_simulation_variables();
    
    // return immediately if no nonlinear terms
    if (simVars->pde.use_linear_div == 1)
      {
	phi_F1.physical_set_zero();
	vort_F1.physical_set_zero();
	div_F1.physical_set_zero();

	return;
      }


    if (simVars->libpfasst.implicit_coriolis_force) 
      {
   
	SWE_Sphere_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());
		  
	// compute the explicit nonlinear right-hand side
	timestepper->euler_timestep_update_nonlinear(
						     phi_Y, 
						     vort_Y,
						     div_Y,
						     phi_F1,
						     vort_F1,
						     div_F1, 
						     simVars->timecontrol.max_simulation_time
						     );
      }   
    else
      {

	if (!simVars->sim.use_topography) 
	  {

	    SWE_Sphere_TS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(i_Y->get_level());
	    
	    // compute the explicit nonlinear right-hand side
	    timestepper->euler_timestep_update_coriolis_and_nonlinear(
								      phi_Y, 
								      vort_Y,
								      div_Y,
								      phi_F1,
								      vort_F1,
								      div_F1, 
								      simVars->timecontrol.max_simulation_time
								      );
	  } 
	else
	  {

	    SWE_Sphere_TS_lg_erk_lc_n_t_erk* timestepper = i_ctx->get_lg_erk_lc_n_t_erk_timestepper(i_Y->get_level());
	    
	    // compute the explicit nonlinear right-hand side
	    timestepper->euler_timestep_update_coriolis_and_nonlinear(
								      phi_Y, 
								      vort_Y,
								      div_Y,
								      phi_F1,
								      vort_F1,
								      div_F1, 
								      simVars->timecontrol.max_simulation_time
								      );

	    
	  }

      }
  }

  // evaluates the first implicit piece o_F2 = F2(i_Y)
  void ceval_f2 (
		 SphereDataVars *i_Y, 
		 double i_t, 
		 SphereDataCtx *i_ctx, 
		 SphereDataVars *o_F2
		 ) 
  {
    const SphereData& phi_Y  = i_Y->get_phi();
    const SphereData& vort_Y = i_Y->get_vort();
    const SphereData& div_Y  = i_Y->get_div();

    SphereData& phi_F2  = o_F2->get_phi();
    SphereData& vort_F2 = o_F2->get_vort();
    SphereData& div_F2  = o_F2->get_div();

    phi_F2.physical_set_zero();
    vort_F2.physical_set_zero();
    div_F2.physical_set_zero();

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    if (simVars->libpfasst.implicit_coriolis_force) 
      {
    
	// get the explicit timestepper 
	SWE_Sphere_TS_l_erk_n_erk* timestepper = i_ctx->get_l_erk_n_erk_timestepper(i_Y->get_level());
	
	// compute the linear right-hand side
	timestepper->euler_timestep_update_linear(
						  phi_Y, 
						  vort_Y,
						  div_Y,
						  phi_F2,
						  vort_F2,
						  div_F2,
						  simVars->timecontrol.max_simulation_time
						  );
      }
    else
      {

	if (!simVars->sim.use_topography) 
	  {

	    // get the explicit timestepper 
	    SWE_Sphere_TS_lg_erk_lc_n_erk* timestepper = i_ctx->get_lg_erk_lc_n_erk_timestepper(i_Y->get_level());
	    
	    // compute the linear right-hand side
	    timestepper->euler_timestep_update_linear(
						      phi_Y, 
						      vort_Y,
						      div_Y,
						      phi_F2,
						      vort_F2,
						      div_F2,
						      simVars->timecontrol.max_simulation_time
						      );

	  }
	else 
	  {

	    // get the explicit timestepper 
	    SWE_Sphere_TS_lg_erk_lc_n_t_erk* timestepper = i_ctx->get_lg_erk_lc_n_t_erk_timestepper(i_Y->get_level());
	    
	    // compute the linear right-hand side
	    timestepper->euler_timestep_update_linear(
						      phi_Y, 
						      vort_Y,
						      div_Y,
						      phi_F2,
						      vort_F2,
						      div_F2,
						      simVars->timecontrol.max_simulation_time
						      );


	  }
      }
  }

  // solves the first implicit system for io_Y
  // then updates o_F2 with the new value of F2(io_Y)
  void ccomp_f2 (
		 SphereDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 SphereDataVars *i_Rhs, 
		 SphereDataCtx *i_ctx, 
		 SphereDataVars *o_F2
		 ) 
  {
    SphereData& phi_Y  = io_Y->get_phi();
    SphereData& vort_Y = io_Y->get_vort();
    SphereData& div_Y  = io_Y->get_div();

    const SphereData& phi_Rhs  = i_Rhs->get_phi();
    const SphereData& vort_Rhs = i_Rhs->get_vort();
    const SphereData& div_Rhs  = i_Rhs->get_div();

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // first copy the rhs into the solution vector
    // this is needed to call the SWEET function run_ts_timestep
    phi_Y  = phi_Rhs;
    vort_Y = vort_Rhs;
    div_Y  = div_Rhs;
    

    if (simVars->libpfasst.use_rexi) 
      {
	// get the rexi time stepper
	SWE_Sphere_TS_l_rexi* timestepper = i_ctx->get_l_rexi_timestepper(io_Y->get_level());

	// solve the implicit system using the Helmholtz solver
	timestepper->run_timestep(
				  phi_Y,
				  vort_Y,
				  div_Y,
				  i_dt,
				  simVars->timecontrol.max_simulation_time
				  );
	
      }
    else 
      {
	if (simVars->libpfasst.implicit_coriolis_force) 
	  {
	    
	    // get the irk timestepper
	    SWE_Sphere_TS_l_irk* timestepper = i_ctx->get_l_irk_timestepper(io_Y->get_level());
	    
	    // solve the implicit system using the Helmholtz solver
	    timestepper->run_timestep(
				      phi_Y,
				      vort_Y,
				      div_Y,
				      i_dt,
				      simVars->timecontrol.max_simulation_time
				      );

	  }
	else
	  {
	    
	    // get the irk timestepper
	    SWE_Sphere_TS_lg_irk* timestepper = i_ctx->get_lg_irk_timestepper(io_Y->get_level());
	    
	    // solve the implicit system using the Helmholtz solver
	    timestepper->run_timestep(
				      phi_Y,
				      vort_Y,
				      div_Y,
				      i_dt,
				      simVars->timecontrol.max_simulation_time
				      );


	  }
      }

    // now recompute F2 with the new value of Y
    ceval_f2(
      	     io_Y, 
     	     i_t, 
     	     i_ctx, 
     	     o_F2
     	     );
    

    SphereData& phi_F2  = o_F2->get_phi();
    SphereData& vort_F2 = o_F2->get_vort();
    SphereData& div_F2  = o_F2->get_div();

    SphereData phi_F2_new(i_ctx->get_sphere_data_config(io_Y->get_level()));
    SphereData vort_F2_new(i_ctx->get_sphere_data_config(io_Y->get_level()));
    SphereData div_F2_new(i_ctx->get_sphere_data_config(io_Y->get_level()));

    phi_F2_new  = (phi_Y  - phi_Rhs)  / i_dt;
    vort_F2_new = (vort_Y - vort_Rhs) / i_dt;
    div_F2_new  = (div_Y  - div_Rhs)  / i_dt;

    // write_file(*i_ctx, phi_F2,  "prog_phi_F2");
    // write_file(*i_ctx, div_F2,  "prog_div_F2");
    // write_file(*i_ctx, vort_F2, "prog_vort_F2");
    
    // write_file(*i_ctx, phi_F2_new,  "prog_phi_F2_new");
    // write_file(*i_ctx, div_F2_new,  "prog_div_F2_new");
    // write_file(*i_ctx, vort_F2_new, "prog_vort_F2_new");

    // write_file(*i_ctx, (phi_F2  - phi_F2_new)/phi_F2_new,  "prog_phi_F2_diff");
    // write_file(*i_ctx, (div_F2  - div_F2_new)/div_F2_new,  "prog_div_F2_diff");
    // write_file(*i_ctx, (vort_F2 - vort_F2_new)/vort_F2_new, "prog_vort_F2_diff");

    phi_F2  = phi_F2_new;
    vort_F2 = vort_F2_new;
    div_F2  = div_F2_new;

  }

  // evaluates the second implicit piece o_F3 = F3(i_Y)
  // currently contains the implicit artificial viscosity term
  void ceval_f3 (
		 SphereDataVars *i_Y, 
		 double i_t, 
		 int i_level,
		 SphereDataCtx *i_ctx, 
		 SphereDataVars *o_F3
		 ) 
  {
    const SphereData& phi_Y  = i_Y->get_phi();
    const SphereData& vort_Y = i_Y->get_vort();
    const SphereData& div_Y  = i_Y->get_div();

    SphereData& phi_F3  = o_F3->get_phi();
    SphereData& vort_F3 = o_F3->get_vort();
    SphereData& div_F3  = o_F3->get_div();

    // initialize F3 to zero in case no artificial viscosity
    c_sweet_data_setval(o_F3, 0.0);

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // no need to do anything if no artificial viscosity
    if (simVars->sim.viscosity == 0) 
      return;

    // get the parameters used to apply diffusion
    const double r    = simVars->sim.earth_radius;
    const double visc = simVars->sim.viscosity;

    phi_F3 = phi_Y;
    phi_F3.spectral_update_lambda(
				  [&](
				      int n, int m,
				      std::complex<double> &io_data
				      )
				  {
				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
				  }
				  );

    vort_F3 = vort_Y;
    vort_F3.spectral_update_lambda(
				  [&](
				      int n, int m,
				      std::complex<double> &io_data
				      )
				  {
				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
				  }
				  );

    div_F3 = div_Y;
    div_F3.spectral_update_lambda(
				  [&](
				      int n, int m,
				      std::complex<double> &io_data
				      )
				  {
				    io_data *= (-visc*(double)n*((double)n+1.0))/(r*r);
				  }
				  );
  }

  // solves the second implicit system for io_Y
  // then updates o_F3 with the new value o_F3 = F3(io_Y)
  // currently solve the implicit system formed with artificial viscosity term
  void ccomp_f3 (
		 SphereDataVars *io_Y, 
		 double i_t, 
		 double i_dt,
		 int i_level,
		 SphereDataVars *i_Rhs, 
		 SphereDataCtx *i_ctx, 
		 SphereDataVars *o_F3
		 ) 
  {
    SphereData& phi_Y  = io_Y->get_phi();
    SphereData& vort_Y = io_Y->get_vort();
    SphereData& div_Y  = io_Y->get_div();

    SphereData& phi_Rhs  = i_Rhs->get_phi();
    SphereData& vort_Rhs = i_Rhs->get_vort();
    SphereData& div_Rhs  = i_Rhs->get_div();

    // initialize F3 to zero in case no artificial viscosity
    c_sweet_data_setval(o_F3, 0.0);
            
    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // no need to do anything if no artificial viscosity
    if (simVars->sim.viscosity == 0)
      return;

    // get the parameters used to apply diffusion
    const double scalar = simVars->sim.viscosity*i_dt;
    const double r      = simVars->sim.earth_radius;

    // solve (1-dt*visc*diff_op)*rhs = y
    phi_Y  = phi_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r); 
    vort_Y = vort_Rhs.spectral_solve_helmholtz(1.0, -scalar, r); 
    div_Y  = div_Rhs.spectral_solve_helmholtz( 1.0, -scalar, r); 
    
    // now recompute F3 with the new value of Y
    
    // ceval_f3(
    //  	     io_Y, 
    // 	     i_t, 
    // 	     i_level,
    // 	     i_ctx, 
    // 	     o_F3
    // 	     );
    
    SphereData& phi_F3  = o_F3->get_phi();
    SphereData& vort_F3 = o_F3->get_vort();
    SphereData& div_F3  = o_F3->get_div();

    phi_F3  = (phi_Y  - phi_Rhs)  / i_dt;
    vort_F3 = (vort_Y - vort_Rhs) / i_dt;
    div_F3  = (div_Y  - div_Rhs)  / i_dt;
  
  }

}
