#include <iomanip>
#include <math.h>
#include <string>

#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_l_irk_n_erk.hpp"
#include "SWE_Plane_TS_ln_erk.hpp"

#include "ceval.hpp"

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

  // Write the data into the file
  const char* filename_template = simVars->misc.output_file_name_prefix.c_str();
  sprintf(buffer, 
	  filename_template, 
	  i_name, 
	  simVars->timecontrol.current_timestep_size);
  i_planeData.file_physical_saveData_ascii(buffer);
  
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
    PlaneData& h_Y = o_Y->get_h();
    PlaneData& u_Y = o_Y->get_u();
    PlaneData& v_Y = o_Y->get_v();

    // set the time stepping params
    i_ctx->setup_time_steps(
			    i_t,
			    i_dt
			    );

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());
    
    h_Y.physical_set_all(simVars->sim.h0);
    u_Y.physical_set_all(0);
    v_Y.physical_set_all(0);
    
    // initialize height
    h_Y.physical_update_lambda_array_indices(
    					     [&](int i, int j, double &io_data)
    					     {
    					       double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
    					       double y = (((double)j+0.5)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
    					       io_data = SWEPlaneBenchmarks::return_h(*simVars, x, y);
    					     }
    					     );
    
    // initialize velocities
    u_Y.physical_update_lambda_array_indices(
    					     [&](int i, int j, double &io_data)
    					     {
    					       double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
    					       double y = (((double)j+0.5)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
    					       io_data = SWEPlaneBenchmarks::return_u(*simVars, x, y);
    					     }
    					     );
    v_Y.physical_update_lambda_array_indices(
    					     [&](int i, int j, double &io_data)
    					     {
    					       double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
    					       double y = (((double)j+0.5)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
    					       io_data = SWEPlaneBenchmarks::return_v(*simVars, x, y);
    					     }
    					     );
  }

  // finalizes the time step when libpfasst is done 
  // currently does nothing else than outputting the solution
  void cfinal(
	      PlaneDataCtx *i_ctx, 
	      PlaneDataVars *i_Y,
	      int i_nnodes, 
	      int i_niters
	      ) 
  {
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());
    
    /*
    std::vector<int> coarsening_factors(5,0);
    coarsening_factors[0] = 2;
    coarsening_factors[1] = 4;
    coarsening_factors[2] = 8;
    coarsening_factors[3] = 16;
    coarsening_factors[4] = 32;
    
    for (unsigned int i = 0; i < coarsening_factors.size(); ++i) {

      PlaneDataConfig coarsened_dataConfig;
      coarsened_dataConfig.setupAutoPhysicalSpace(
						  simVars->disc.res_spectral[0]/coarsening_factors[i],
						  simVars->disc.res_spectral[1]/coarsening_factors[i]
						  );
      coarsened_dataConfig.printInformation();
      PlaneData coarsened_h_Y 
	= h_Y.spectral_returnWithDifferentModes(&coarsened_dataConfig);

      std::string filename = "prog_h_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters)+"_coarsening_factor_"+std::to_string(coarsening_factors[i]);
      write_file(*i_ctx, coarsened_h_Y, filename.c_str());
      
      PlaneData coarsened_u_Y 
	= u_Y.spectral_returnWithDifferentModes(&coarsened_dataConfig);

      filename = "prog_u_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters)+"_coarsening_factor_"+std::to_string(coarsening_factors[i]);
      write_file(*i_ctx, coarsened_u_Y, filename.c_str());

      PlaneData coarsened_v_Y 
	= v_Y.spectral_returnWithDifferentModes(&coarsened_dataConfig);
      
      filename = "prog_v_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters)+"_coarsening_factor_"+std::to_string(coarsening_factors[i]);
      write_file(*i_ctx, coarsened_v_Y, filename.c_str());
      
    }
    */
    
    std::string filename = "prog_h_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, h_Y, filename.c_str());
    
    filename = "prog_u_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, u_Y, filename.c_str());

    filename = "prog_v_nnodes_"+std::to_string(i_nnodes)+"_niters_"+std::to_string(i_niters);
    write_file(*i_ctx, v_Y, filename.c_str());
    
  }

  // computes a reference solution to check libpfasst's results
  // based on the fourth order explicit (linear and nonlinear) RK
  void creference( 
		  double i_t,
		  PlaneDataCtx *i_ctx,
		  PlaneDataVars *o_Y
		   )
  {
    // get the time step parameters
    SimulationVariables* simVars = i_ctx->get_simulation_variables();
    
    // compute the initial condition
    cinitial(
	     i_ctx, 
	     simVars->timecontrol.current_simulation_time,
	     simVars->timecontrol.current_timestep_size,
	     o_Y
	     );
    
    PlaneData& h = o_Y->get_h();
    PlaneData& u = o_Y->get_u();
    PlaneData& v = o_Y->get_v();
    
    // get the reference timestepper
    SWE_Plane_TS_ln_erk* timestepper = i_ctx->get_reference_timestepper(o_Y->get_level());
    
    const int dt_factor = 1; 
    double dt           = simVars->timecontrol.current_timestep_size 
                        / (double)dt_factor;

    double current_simulation_time = 0;
    int nsteps = 0;

    // compute the reference solution (i.e., obtained with the reference time stepper)
    while (current_simulation_time < i_t)
      {
	if (nsteps%100 == 0)
	  std::cout << "current_simulation_time = "
		    << current_simulation_time
		    << std::endl;
  
	timestepper->run_timestep(
				  h,
				  u,
				  v,
				  dt,
				  dt,
				  current_simulation_time,
				  i_t
				  );
	current_simulation_time += dt;
	nsteps += 1;
      }

    write_file(*i_ctx, h, "prog_h_ref");
    write_file(*i_ctx, u, "prog_u_ref");
    write_file(*i_ctx, v, "prog_v_ref");

  }


  // evaluates the explicit (nonlinear) piece
  void ceval_f1(PlaneDataVars *i_Y,
		double i_t,
		PlaneDataCtx *i_ctx,
		PlaneDataVars *o_F1
		)
  {       
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();

    PlaneData& h_F1 = o_F1->get_h();
    PlaneData& u_F1 = o_F1->get_u();
    PlaneData& v_F1 = o_F1->get_v();
    
    // get the time step parameters
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // get the timestepper
    SWE_Plane_TS_l_irk_n_erk* timestepper = i_ctx->get_l_irk_n_erk_timestepper(i_Y->get_level());
		  
    // compute the explicit nonlinear right-hand side
    timestepper->euler_timestep_update_nonlinear(
						 h_Y, 
						 u_Y,
						 v_Y,
						 h_F1,
						 u_F1,
						 v_F1
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
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();

    PlaneData& h_F2 = o_F2->get_h();
    PlaneData& u_F2 = o_F2->get_u();
    PlaneData& v_F2 = o_F2->get_v();

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();
    
    // get the implicit timestepper 
    SWE_Plane_TS_l_erk* timestepper    = i_ctx->get_l_erk_timestepper(i_Y->get_level());

    // compute the linear right-hand side
    timestepper->euler_timestep_update(
				       h_Y, 
				       u_Y,
				       v_Y,
				       h_F2,
				       u_F2,
				       v_F2,
				       simVars->timecontrol.current_timestep_size,
				       simVars->timecontrol.current_timestep_size,
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
    PlaneData& h_Y = io_Y->get_h();
    PlaneData& u_Y = io_Y->get_u();
    PlaneData& v_Y = io_Y->get_v();

    const PlaneData& h_Rhs = i_Rhs->get_h();
    const PlaneData& u_Rhs = i_Rhs->get_u();
    const PlaneData& v_Rhs = i_Rhs->get_v();

    // first copy the rhs into the solution vector
    // this is needed to call the SWEET function run_ts_timestep
    h_Y = h_Rhs;
    u_Y = u_Rhs;
    v_Y = v_Rhs;

    // get the implicit timestepper 
    SWE_Plane_TS_l_irk_n_erk* timestepper    = i_ctx->get_l_irk_n_erk_timestepper(io_Y->get_level());
    SWE_Plane_TS_l_irk& implicit_timestepper = timestepper->get_implicit_timestepper();

    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    // solve the implicit system using the Helmholtz solver
    implicit_timestepper.run_timestep(
				      h_Y,
				      u_Y,
				      v_Y,
				      i_dt,
				      i_dt,
				      simVars->timecontrol.current_simulation_time,
				      simVars->timecontrol.max_simulation_time
				      );
    
    // now recompute F2 with the new value of Y
    ceval_f2(
	     io_Y, 
	     i_t, 
	     i_ctx, 
	     o_F2
	     );
  }

  // evaluates the second implicit piece o_F3 = F3(i_Y)
  void ceval_f3 (
		 PlaneDataVars *i_Y, 
		 double i_t, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    // not implemented 
  }

  // solves the second implicit system for io_Y
  // then updates o_F3 with the new value o_F3 = F3(io_Y)
  void ccomp_f3 (
		 PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 PlaneDataVars *i_Rhs, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    // not implemented
  }

  // applies implicit viscosity to the variables
  void capply_viscosity(PlaneDataVars *io_Y,
			double i_t, 
			double i_dt,
			int i_level, 
			PlaneDataCtx *i_ctx
			)
  {
    // get the simulation variables
    SimulationVariables* simVars = i_ctx->get_simulation_variables();

    if (simVars->sim.viscosity == 0)
      return;
		  
    PlaneData& h_Y = io_Y->get_h();
    PlaneData& u_Y = io_Y->get_u();
    PlaneData& v_Y = io_Y->get_v();

    // get the operators for this level
    PlaneOperators* op = i_ctx->get_plane_operators(i_level);

    // use the implicit diffusion from SWEET
    h_Y = op->implicit_diffusion(h_Y, i_dt*simVars->sim.viscosity, simVars->sim.viscosity_order );
    u_Y = op->implicit_diffusion(u_Y, i_dt*simVars->sim.viscosity, simVars->sim.viscosity_order );
    v_Y = op->implicit_diffusion(v_Y, i_dt*simVars->sim.viscosity, simVars->sim.viscosity_order );

  }
}
