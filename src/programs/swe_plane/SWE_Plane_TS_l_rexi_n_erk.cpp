/*
 * SWE_Plane_TS_l_rexi_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "../swe_plane/SWE_Plane_TS_l_rexi_n_erk.hpp"

#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/SimulationVariables.hpp>



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_rexi_n_erk::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double i_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u*op.diff_c_x(i_h) + i_v*op.diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);

}


void SWE_Plane_TS_l_rexi_n_erk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_rexi_n_erk: Only constant time step size allowed");

	if (timestepping_order_nonlinear == 1)
	{
		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				i_dt,
				i_simulation_timestamp
			);

		// standard time stepping
		timestepping_rk.run_timestep(
				this,
				&SWE_Plane_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				timestepping_order_nonlinear,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order_nonlinear == 2)
	{
		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				i_dt*0.5,
				i_simulation_timestamp
			);

		// standard time stepping
		timestepping_rk.run_timestep(
				this,
				&SWE_Plane_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				timestepping_order_nonlinear,
				i_simulation_timestamp
			);

		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				i_dt*0.5,
				i_simulation_timestamp
			);
	}
	else
	{
		FatalError("SWE_Plane_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_n_erk::setup(
		REXI_SimulationVariables &i_rexi,

		int i_nonlinear_order,

		bool i_use_only_linear_divergence
)
{
	use_only_linear_divergence = i_use_only_linear_divergence;

	ts_l_rexi.setup(i_rexi, "phi0", simVars.timecontrol.current_timestep_size);

	timestepping_order_nonlinear = i_nonlinear_order;
	timestepping_rk.setupBuffers(op.planeDataConfig, timestepping_order_nonlinear);

	if (simVars.disc.space_grid_use_c_staggering)
		FatalError("Staggering not supported for l_rexi_n_erk");

}


SWE_Plane_TS_l_rexi_n_erk::SWE_Plane_TS_l_rexi_n_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_l_rexi(simVars, op)
{
}



SWE_Plane_TS_l_rexi_n_erk::~SWE_Plane_TS_l_rexi_n_erk()
{
}

