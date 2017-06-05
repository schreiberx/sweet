/*
 * SWE_Plane_TS_l_rexi_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_rexi_n_erk.hpp"
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>



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

}


void SWE_Plane_TS_l_rexi_n_erk::run_timestep(
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
		FatalError("SWE_Plane_TS_l_rexi_n_erk: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;

	o_dt = i_fixed_dt;

	ts_l_rexi.run_timestep(
			io_h, io_u, io_v,
			o_dt,
			i_fixed_dt,
			i_simulation_timestamp,
			i_max_simulation_time
		);

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Plane_TS_l_rexi_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			o_dt,
			i_fixed_dt,
			timestepping_order_nonlinear,
			i_simulation_timestamp,
			i_max_simulation_time
		);
}



/*
 * Setup
 */
void SWE_Plane_TS_l_rexi_n_erk::setup(
		double i_h,						///< sampling size
		int i_M,						///< number of sampling points
		int i_L,						///< number of sampling points for Gaussian approximation
										///< set to 0 for auto detection

		bool i_rexi_half,				///< use half-pole reduction
		bool i_rexi_normalization,		///< REXI normalization

		int i_nonlinear_order
)
{
	ts_l_rexi.setup(i_h, i_M, i_L, i_rexi_half, i_rexi_normalization);

	timestepping_order_nonlinear = i_nonlinear_order;
	timestepping_rk.setupBuffers(op.planeDataConfig, timestepping_order_nonlinear);
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

