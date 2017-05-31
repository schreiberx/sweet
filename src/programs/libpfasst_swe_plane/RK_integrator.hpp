#ifndef _RK_INTEGRATOR_HPP_
#define _RK_INTEGRATOR_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>

#include <sweet/SimulationVariables.hpp>

#include "LevelSingleton.hpp"

// Helper class to launch the RK integrator implemented in sweet
// this is used to compute a reference solution to check Libpfasst's results

class RK_integrator
{

public:
  
  // constructor
  RK_integrator(
                SimulationVariables *i_simVars,
                PlaneOperators *i_op
	       )
    : simVars(i_simVars),
      op(i_op)
  {}

  /**
   * Execute a single simulation time step
   */
  void p_run_timestep(
		      PlaneData &io_prog_phi,
		      PlaneData &io_prog_u,
		      PlaneData &io_prog_v,
		      double i_dt
		      )
  {
    // select fourth order runge-kutta integrator
    simVars->disc.timestepping_order = 4; 
    // timestepping parameters
    simVars->timecontrol.current_timestep_size   = i_dt;
    simVars->timecontrol.current_simulation_time = 0.0;
    simVars->timecontrol.max_simulation_time     = i_dt;

    // standard time stepping
    timestepping.run_timestep(
			      this,
			      &RK_integrator::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
			      io_prog_phi, io_prog_u, io_prog_v,
			      i_dt,
			      i_dt,
			      simVars->disc.timestepping_order,
			      simVars->timecontrol.current_simulation_time,
			      simVars->timecontrol.max_simulation_time
			      );
  }
    
private:
  
  // timestepping object containing the parameters for RK
  PlaneDataTimesteppingRK timestepping;
  // pointer to the SimulationVariables oject
  SimulationVariables    *simVars;
  // pointer to the LevelSingleton object for this level
  PlaneOperators         *op;
  
  /**
   * Main routine for method to be used in case of finite differences
   */
  
  void p_run_euler_timestep_update(
				   const PlaneData &i_phi,	///< prognostic variables
				   const PlaneData &i_u,	///< prognostic variables
				   const PlaneData &i_v,	///< prognostic variables
				   
				   PlaneData &o_phi_t,	///< time updates
				   PlaneData &o_u_t,	///< time updates
				   PlaneData &o_v_t,	///< time updates
				   
				   double &o_dt,			///< time step restriction
				   double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
				   double i_simulation_timestamp = -1
				   )
  {
    /**
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * TODO: rearrange equations to phi formulation
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     */
    o_dt = simVars->timecontrol.current_timestep_size;
    
    /*
     * linearized non-conservative (advective) formulation:
     *
     * phi_t = -(h0*g)*u_x - (h0*g)*v_ym
     * u_t = -phi_x + f*v
     * v_t = -phi_y - f*u
     */
    
    o_u_t = -op->diff_c_x(i_phi);
    o_v_t = -op->diff_c_y(i_phi);
    
    o_u_t += simVars->sim.f0*i_v;
    o_v_t -= simVars->sim.f0*i_u;
    
    // standard update
    o_phi_t = -(op->diff_c_x(i_u) + op->diff_c_y(i_v))*(simVars->sim.h0*simVars->sim.gravitation);
    
    assert(simVars->sim.viscosity_order == 2);
    if (simVars->sim.viscosity != 0)
      {
	o_phi_t += op->laplace(i_phi);
	o_u_t += op->laplace(i_u);
	o_v_t += op->laplace(i_v);
      }
    
  }
  
  // default constructor, copy constructor, and operator= are disabled
  RK_integrator() {};
  RK_integrator(const RK_integrator&);
  RK_integrator& operator=(const RK_integrator&);
};

#endif
