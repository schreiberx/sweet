#include <iomanip>
#include <math.h>

#include <sweet/SimulationVariables.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>

#include "ceval.hpp"
#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_l_direct.hpp"

const double adv_coef  = 0.0;
const double diff_coef = 0.0;
const double reac_coef = 0.0;

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
    
    // initialize geopotential
    h_Y.physical_update_lambda_array_indices(
					     [&](int i, int j, double &io_data)
					     {
					       double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
					       double y = (((double)j+0.5)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
					       
					       io_data = SWEPlaneBenchmarks::return_h(*simVars, x, y);
					     }
					     );
    
    // initialize velocity
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
    
    h_Y.print_physicalArrayData();
    u_Y.print_physicalArrayData();
    v_Y.print_physicalArrayData();

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
    
    std::cerr <<  std::endl;
    h_Y.print_physicalArrayData();
    std::cerr <<  std::endl;
    u_Y.print_physicalArrayData();
    std::cerr <<  std::endl;
    v_Y.print_physicalArrayData();


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
    
    // get the reference timestepper (will be improved later)
    SWE_Plane_TS_interface* timestepper         = i_ctx->get_reference_timestepper(o_Y->get_level());
    SWE_Plane_TS_l_direct* l_direct_timestepper = dynamic_cast<SWE_Plane_TS_l_direct*>(timestepper); 

    // compute the reference solution (i.e., obtained with the reference time stepper)
    l_direct_timestepper->run_timestep(
				       h,
				       u,
				       v,
				       simVars->timecontrol.current_timestep_size,
				       simVars->timecontrol.current_timestep_size,
				       simVars->timecontrol.current_simulation_time,
				       simVars->timecontrol.max_simulation_time
				       );
    
    std::cerr <<  std::endl;
    h.print_physicalArrayData();
    std::cerr <<  std::endl;
    u.print_physicalArrayData();
    std::cerr <<  std::endl;
    v.print_physicalArrayData();
    
  }


  // evaluates the explicit piece
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

    // get the timestepper (will be improved later)
    SWE_Plane_TS_interface* timestepper    = i_ctx->get_timestepper(i_Y->get_level());
    SWE_Plane_TS_l_erk* l_erk_timestepper = dynamic_cast<SWE_Plane_TS_l_erk*>(timestepper); 
		  
    // compute the explicit right-hand side
    l_erk_timestepper->euler_timestep_update(
					     h_Y, 
					     u_Y,
					     v_Y,
					     h_F1,
					     u_F1,
					     v_F1,
					     simVars->timecontrol.current_timestep_size,
					     simVars->timecontrol.current_timestep_size,
					     simVars->timecontrol.max_simulation_time
					     );

  }
  
  // evaluates the first implicit piece
  // currently does nothing since diff_coef = 0
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

    h_F2 = diff_coef * h_Y;
    u_F2 = diff_coef * u_Y;
    v_F2 = diff_coef * v_Y;
  }

  // solves the first implicit system
  // currently does nothing since diff_coef = 0
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

    PlaneData& h_F2 = o_F2->get_h();
    PlaneData& u_F2 = o_F2->get_u();
    PlaneData& v_F2 = o_F2->get_v();

    h_Y = h_Rhs / (1.0 - diff_coef * i_dt);
    u_Y = u_Rhs   / (1.0 - diff_coef * i_dt);
    v_Y = v_Rhs   / (1.0 - diff_coef * i_dt);

    h_F2 = (h_Y - h_Rhs) / i_dt;
    u_F2 = (u_Y   - u_Rhs)   / i_dt;
    v_F2 = (v_Y   - v_Rhs)   / i_dt;

  }

  // evaluates the second implicit piece
  // currently does nothing since reac_coef = 0
  void ceval_f3 (
		 PlaneDataVars *i_Y, 
		 double i_t, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();

    PlaneData& h_F3 = o_F3->get_h();
    PlaneData& u_F3 = o_F3->get_u();
    PlaneData& v_F3 = o_F3->get_v();

    h_F3 = reac_coef * h_Y;
    u_F3 = reac_coef * u_Y;
    v_F3 = reac_coef * v_Y;
  }

  // solves the second implicit system
  // currently does nothing since diff_coef = 0
  void ccomp_f3 (
		 PlaneDataVars *io_Y, 
		 double i_t, 
		 double i_dt, 
		 PlaneDataVars *i_Rhs, 
		 PlaneDataCtx *i_ctx, 
		 PlaneDataVars *o_F3
		 ) 
  {
    PlaneData& h_Y = io_Y->get_h();
    PlaneData& u_Y = io_Y->get_u();
    PlaneData& v_Y = io_Y->get_v();

    const PlaneData& h_Rhs = i_Rhs->get_h();
    const PlaneData& u_Rhs = i_Rhs->get_u();
    const PlaneData& v_Rhs = i_Rhs->get_v();

    PlaneData& h_F3 = o_F3->get_h();
    PlaneData& u_F3 = o_F3->get_u();
    PlaneData& v_F3 = o_F3->get_v();

    h_Y = h_Rhs / (1.0 - reac_coef * i_dt);
    u_Y = u_Rhs   / (1.0 - reac_coef * i_dt);
    v_Y = v_Rhs   / (1.0 - reac_coef * i_dt);

    h_F3 = (h_Y - h_Rhs) / i_dt;
    u_F3 = (u_Y   - u_Rhs)   / i_dt;
    v_F3 = (v_Y   - v_Rhs)   / i_dt;
  }
}
