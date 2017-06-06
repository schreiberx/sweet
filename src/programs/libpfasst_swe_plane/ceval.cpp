#include <iomanip>
#include <math.h>

#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_l_irk_n_erk.hpp"
#include "SWE_Plane_TS_l_direct.hpp"
#include "SWE_Plane_TS_l_rexi.hpp"
#include "SWE_Plane_TS_l_rexi_n_erk.hpp"

#include "ceval.hpp"

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
    
    // Gaussian dam break
    simVars->setup.benchmark_scenario_id = 1;

    h_Y.physical_set_all(simVars->sim.h0);
    u_Y.physical_set_all(0);
    v_Y.physical_set_all(0);
    
    // initialize height
    h_Y.physical_update_lambda_array_indices(
					     [&](int i, int j, double &io_data)
					     {
					       double x = (((double)i)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
					       double y = (((double)j)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
					       io_data = SWEPlaneBenchmarks::return_h(*simVars, x, y);
					     }
					     );
    
    // initialize velocity
    u_Y.physical_update_lambda_array_indices(
					     [&](int i, int j, double &io_data)
					     {
					       double x = (((double)i)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
					       double y = (((double)j)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
					       io_data = SWEPlaneBenchmarks::return_u(*simVars, x, y);
					     }
					     );
    v_Y.physical_update_lambda_array_indices(
					     [&](int i, int j, double &io_data)
					     {
					       double x = (((double)i)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
					       double y = (((double)j)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
					       io_data = SWEPlaneBenchmarks::return_v(*simVars, x, y);
					     }
					     );
    /*
    h_Y.print_physicalArrayData();
    u_Y.print_physicalArrayData();
    v_Y.print_physicalArrayData();
    */
  }

  // finalizes the time step when libpfasst is done 
  // currently does nothing else than outputting the solution
  void cfinal(
	      PlaneDataCtx *i_ctx, 
	      PlaneDataVars *i_Y
	      ) 
  {
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();
    /*
    std::cerr <<  std::endl;
    h_Y.print_physicalArrayData();
    std::cerr <<  std::endl;
    u_Y.print_physicalArrayData();
    std::cerr <<  std::endl;
    v_Y.print_physicalArrayData();
    */

    std::cout << "cfinal is not implemented yet" << std::endl;
    //not implemented yet
  }

  // computes a reference solution to check libpfasst's results
  // based on the "direct solution"
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
    SWE_Plane_TS_l_rexi_n_erk* timestepper = i_ctx->get_reference_timestepper(o_Y->get_level());

    // compute the reference solution (i.e., obtained with the reference time stepper)
    timestepper->run_timestep(
			      h,
			      u,
			      v,
			      simVars->timecontrol.current_timestep_size,
			      simVars->timecontrol.current_timestep_size,
			      simVars->timecontrol.current_simulation_time,
			      simVars->timecontrol.max_simulation_time
			      );
    
    /*std::cerr <<  std::endl;
    h.print_physicalArrayData();
    std::cerr <<  std::endl;
    u.print_physicalArrayData();
    std::cerr <<  std::endl;
    v.print_physicalArrayData();*/
    
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
    SWE_Plane_TS_l_irk_n_erk* timestepper = i_ctx->get_timestepper(i_Y->get_level());
		  
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
    
    // get the simulation parameters
    const double eta_bar = simVars->sim.h0;
    const double g       = simVars->sim.gravitation;
    const double f0      = simVars->sim.f0;

    // get the implicit timestepper 
    SWE_Plane_TS_l_irk_n_erk* timestepper    = i_ctx->get_timestepper(i_Y->get_level());
    SWE_Plane_TS_l_irk& implicit_timestepper = timestepper->get_implicit_timestepper();

    // get the spectral space operators
    PlaneOperatorsComplex& opComplex = implicit_timestepper.get_plane_operators_complex();

    // convert to spectral data
    const PlaneDataComplex eta0_Y = Convert_PlaneData_To_PlaneDataComplex::physical_convert(h_Y);
    const PlaneDataComplex u0_Y   = Convert_PlaneData_To_PlaneDataComplex::physical_convert(u_Y);
    const PlaneDataComplex v0_Y   = Convert_PlaneData_To_PlaneDataComplex::physical_convert(v_Y);
    PlaneDataComplex eta0_F2      = Convert_PlaneData_To_PlaneDataComplex::physical_convert(h_F2);
    PlaneDataComplex u0_F2        = Convert_PlaneData_To_PlaneDataComplex::physical_convert(u_F2);
    PlaneDataComplex v0_F2        = Convert_PlaneData_To_PlaneDataComplex::physical_convert(v_F2);

    // compute the right-hand side for eta
    eta0_F2 = -eta_bar * (opComplex.diff_c_x(u0_Y) + opComplex.diff_c_y(v0_Y));
    
    // compute the right-hand side for u and v
    u0_F2   = -g * opComplex.diff_c_x(eta0_Y) + f0 * v0_Y;
    v0_F2   = -g * opComplex.diff_c_y(eta0_Y) - f0 * u0_Y;
    
    // convert back to physical data
    h_F2 = Convert_PlaneDataComplex_To_PlaneData::physical_convert(eta0_F2);
    u_F2 = Convert_PlaneDataComplex_To_PlaneData::physical_convert(u0_F2);
    v_F2 = Convert_PlaneDataComplex_To_PlaneData::physical_convert(v0_F2);

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
    SWE_Plane_TS_l_irk_n_erk* timestepper    = i_ctx->get_timestepper(io_Y->get_level());
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
				      simVars->timecontrol.current_simulation_time + i_dt
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
}
