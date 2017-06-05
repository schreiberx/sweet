/*
 * SWE_Plane_TS_l_irk1_n_erk1.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_irk1_n_erk1.hpp"







/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_irk1_n_erk1::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t	///< time updates
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	PlaneData hdiff = i_h-simVars.sim.h0;

	o_h_t = -op.diff_c_x(i_u*hdiff) - op.diff_c_y(i_v*hdiff);
	o_u_t = -simVars.sim.gravitation*op.diff_c_x(hdiff) - i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -simVars.sim.gravitation*op.diff_c_y(hdiff) - i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

}


void SWE_Plane_TS_l_irk1_n_erk1::run_timestep(
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
		FatalError("Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;

	o_dt = i_fixed_dt;


	PlaneData h_linear_t1 = io_h;
	PlaneData u_linear_t1 = io_u;
	PlaneData v_linear_t1 = io_v;

	ts_l_irk.run_timestep(
			h_linear_t1, u_linear_t1, v_linear_t1,
			o_dt,
			i_fixed_dt,
			i_simulation_timestamp,
			i_max_simulation_time
		);

	// compute non-linear tendencies at half time step
	PlaneData h_dt_nonlinear(op.planeDataConfig);
	PlaneData u_dt_nonlinear(op.planeDataConfig);
	PlaneData v_dt_nonlinear(op.planeDataConfig);

#if 1

	// standard time stepping
	euler_timestep_update_nonlinear(
			io_h, io_u, io_v,
			h_dt_nonlinear, u_dt_nonlinear, v_dt_nonlinear
		);

#else
	// compute value at half time step
	PlaneData h_dh_nonlinear = (io_h + h_linear_t1)*0.5;
	PlaneData u_dh_nonlinear = (io_u + u_linear_t1)*0.5;
	PlaneData v_dh_nonlinear = (io_v + v_linear_t1)*0.5;

	// MIDPOINT
	euler_timestep_update_nonlinear(
			h_dh_nonlinear, u_dh_nonlinear, v_dh_nonlinear,
			h_dt_nonlinear, u_dt_nonlinear, v_dt_nonlinear
		);

#endif


	io_h = h_linear_t1 + h_dt_nonlinear*i_fixed_dt;
	io_u = u_linear_t1 + u_dt_nonlinear*i_fixed_dt;
	io_v = v_linear_t1 + v_dt_nonlinear*i_fixed_dt;
}



/*
 * Setup
 */
void SWE_Plane_TS_l_irk1_n_erk1::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
	ts_l_irk.setup(timestepping_order);

	if (timestepping_order != 1)
		FatalError("SWE_Plane_TS_l_irk1_n_erk1: Only 1st order TS supported with this implementation");
}


SWE_Plane_TS_l_irk1_n_erk1::SWE_Plane_TS_l_irk1_n_erk1(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_l_irk(simVars, op)
{
	setup(simVars.disc.timestepping_order);
	ts_l_irk.setup(simVars.disc.timestepping_order);
}



SWE_Plane_TS_l_irk1_n_erk1::~SWE_Plane_TS_l_irk1_n_erk1()
{
}

