/*
 * SWE_Plane_TS_l_phi0_n_edt.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_phi0_n_edt.hpp"



/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_phi0_n_edt::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double &o_dt,
		double i_fixed_dt,
		double i_max_timestamp
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

	o_dt = i_fixed_dt;
}



void SWE_Plane_TS_l_phi0_n_edt::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,
		double i_max_simulation_time
)
{
	if (i_fixed_dt <= 0)
		FatalError("SWE_Plane_TS_l_phi0_n_edt: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;

	int timestepping_order_nonlinear = 1;
	if (timestepping_order_nonlinear == 1)
	{
		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				o_dt,
				i_fixed_dt,
				i_simulation_timestamp,
				i_max_simulation_time
			);
	}
	else if (timestepping_order_nonlinear == 2)
	{
		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				o_dt,
				i_fixed_dt*0.5,
				i_simulation_timestamp,
				i_max_simulation_time
			);
#if 0
		// standard time stepping
		timestepping_rk.run_timestep(
				this,
				&SWE_Plane_TS_l_phi0_n_edt::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				o_dt,
				i_fixed_dt,
				timestepping_order_nonlinear,
				i_simulation_timestamp,
				i_max_simulation_time
			);
#endif
		ts_l_rexi.run_timestep(
				io_h, io_u, io_v,
				o_dt,
				i_fixed_dt*0.5,
				i_simulation_timestamp,
				i_max_simulation_time
			);
	}

	o_dt = i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_phi0_n_edt::setup(
		REXI_SimulationVariables &i_rexiSimVars
)
{
	ts_l_rexi.setup(i_rexiSimVars);
}


SWE_Plane_TS_l_phi0_n_edt::SWE_Plane_TS_l_phi0_n_edt(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_l_rexi(simVars, op)
{
	FatalError("lkjasdf");
}



SWE_Plane_TS_l_phi0_n_edt::~SWE_Plane_TS_l_phi0_n_edt()
{
}

