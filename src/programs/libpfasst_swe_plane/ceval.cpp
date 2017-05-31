#include <iomanip>
#include <math.h>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>
#include "ceval.hpp"
#include "RK_integrator.hpp"

const double adv_coef  = 0.0;
const double diff_coef = 0.0;
const double reac_coef = 0.0;

extern "C"
{
  // initialization of the variables (initial condition)
  void cinitial(
		PlaneDataCtx *i_ctx, 
		PlaneDataVars *o_Y
		) 
  {
    PlaneData& phi_Y = o_Y->get_phi();
    PlaneData& u_Y   = o_Y->get_u();
    PlaneData& v_Y   = o_Y->get_v();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());
    
    // Gaussian dam break
    simVars->setup.benchmark_scenario_id = 1;

    phi_Y.physical_set_all(simVars->sim.h0 * simVars->sim.gravitation);
    u_Y.physical_set_all(0);
    v_Y.physical_set_all(0);
    
    // initialize geopotential
    phi_Y.physical_update_lambda_array_indices(
					       [&](int i, int j, double &io_data)
					       {
						 double x = (((double)i+0.5)/(double)simVars->disc.res_physical[0])*simVars->sim.domain_size[0];
						 double y = (((double)j+0.5)/(double)simVars->disc.res_physical[1])*simVars->sim.domain_size[1];
						 
						 io_data = SWEPlaneBenchmarks::return_h(*simVars, x, y);
					       }
					       );
    phi_Y *= simVars->sim.gravitation;

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
    
    phi_Y.print_physicalArrayData();
    u_Y.print_physicalArrayData();
    v_Y.print_physicalArrayData();

  }

  // finalizes the time step when libpfasst is done 
  // currently does nothing
  void cfinal(
	      PlaneDataCtx *i_ctx, 
	      PlaneDataVars *i_Y
	      ) 
  {
    std::cout << "cfinal is not implemented yet" << std::endl;
    //not implemented yet
  }

  // computes a reference solution to check libpfasst's results
  // currently based on another sweet integrator (RK) 
  // but later based on the "direct solution"
  void creference( 
	      double i_t,
	      PlaneDataCtx *i_ctx,
	      PlaneDataVars *o_Y
	      )
  {
    // compute the initial condition
    cinitial(
	     i_ctx, 
	     o_Y
	     );
    
    PlaneData& phi = o_Y->get_phi();
    PlaneData& u   = o_Y->get_u();
    PlaneData& v   = o_Y->get_v();
    
    // instante and setup the RK integrator from sweet
    RK_integrator rki(
		      i_ctx->get_simulation_variables(),
		      i_ctx->get_plane_operators(o_Y->get_level())
		      );

    // compute the solution with RK
    rki.p_run_timestep(
		       phi,
		       u,
		       v,
		       i_t
		       );
    
    std::cerr <<  std::endl;
    phi.print_physicalArrayData();
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
    const PlaneData& phi_Y = i_Y->get_phi();
    const PlaneData& u_Y   = i_Y->get_u();
    const PlaneData& v_Y   = i_Y->get_v();

    PlaneData& phi_F1 = o_F1->get_phi();
    PlaneData& u_F1   = o_F1->get_u();
    PlaneData& v_F1   = o_F1->get_v();

    /*
     * linearized non-conservative (advective) formulation:
     *
     * phi_t = -(h0*g)*u_x - (h0*g)*v_ym
     * u_t = -phi_x + f*v
     * v_t = -phi_y - f*u
     */
    
    SimulationVariables *simVars = i_ctx->get_simulation_variables();
    PlaneOperators      *op      = i_ctx->get_plane_operators(i_Y->get_level());

    u_F1 = -op->diff_c_x(phi_Y);
    v_F1 = -op->diff_c_y(phi_Y);

    u_F1 += simVars->sim.f0 * v_Y;
    v_F1 -= simVars->sim.f0 * u_Y;

    // standard update
    phi_F1 = -(op->diff_c_x(u_Y) + op->diff_c_y(v_Y))
           * simVars->sim.h0 * simVars->sim.gravitation;

    assert(simVars->sim.viscosity_order == 2);
    if (simVars->sim.viscosity != 0)
      {
    	phi_F1 += op->laplace(phi_Y);
    	u_F1   += op->laplace(u_Y);
    	v_F1   += op->laplace(v_Y);
      }
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
    const PlaneData& phi_Y = i_Y->get_phi();
    const PlaneData& u_Y   = i_Y->get_u();
    const PlaneData& v_Y   = i_Y->get_v();

    PlaneData& phi_F2 = o_F2->get_phi();
    PlaneData& u_F2   = o_F2->get_u();
    PlaneData& v_F2   = o_F2->get_v();

    phi_F2 = diff_coef * phi_Y;
    u_F2   = diff_coef * u_Y;
    v_F2   = diff_coef * v_Y;
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
    PlaneData& phi_Y = io_Y->get_phi();
    PlaneData& u_Y   = io_Y->get_u();
    PlaneData& v_Y   = io_Y->get_v();

    const PlaneData& phi_Rhs = i_Rhs->get_phi();
    const PlaneData& u_Rhs   = i_Rhs->get_u();
    const PlaneData& v_Rhs   = i_Rhs->get_v();

    PlaneData& phi_F2 = o_F2->get_phi();
    PlaneData& u_F2   = o_F2->get_u();
    PlaneData& v_F2   = o_F2->get_v();

    phi_Y = phi_Rhs / (1.0 - diff_coef * i_dt);
    u_Y   = u_Rhs   / (1.0 - diff_coef * i_dt);
    v_Y   = v_Rhs   / (1.0 - diff_coef * i_dt);

    phi_F2 = (phi_Y - phi_Rhs) / i_dt;
    u_F2   = (u_Y   - u_Rhs)   / i_dt;
    v_F2   = (v_Y   - v_Rhs)   / i_dt;

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
    const PlaneData& phi_Y = i_Y->get_phi();
    const PlaneData& u_Y   = i_Y->get_u();
    const PlaneData& v_Y   = i_Y->get_v();

    PlaneData& phi_F3 = o_F3->get_phi();
    PlaneData& u_F3   = o_F3->get_u();
    PlaneData& v_F3   = o_F3->get_v();

    phi_F3 = reac_coef * phi_Y;
    u_F3   = reac_coef * u_Y;
    v_F3   = reac_coef * v_Y;
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
    PlaneData& phi_Y = io_Y->get_phi();
    PlaneData& u_Y   = io_Y->get_u();
    PlaneData& v_Y   = io_Y->get_v();

    const PlaneData& phi_Rhs = i_Rhs->get_phi();
    const PlaneData& u_Rhs   = i_Rhs->get_u();
    const PlaneData& v_Rhs   = i_Rhs->get_v();

    PlaneData& phi_F3 = o_F3->get_phi();
    PlaneData& u_F3   = o_F3->get_u();
    PlaneData& v_F3   = o_F3->get_v();

    phi_Y = phi_Rhs / (1.0 - reac_coef * i_dt);
    u_Y   = u_Rhs   / (1.0 - reac_coef * i_dt);
    v_Y   = v_Rhs   / (1.0 - reac_coef * i_dt);

    phi_F3 = (phi_Y - phi_Rhs) / i_dt;
    u_F3   = (u_Y   - u_Rhs)   / i_dt;
    v_F3   = (v_Y   - v_Rhs)   / i_dt;
  }
}
