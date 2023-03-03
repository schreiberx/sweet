/*
 * SWE_Plane_TS_l_irk_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */


#include "SWE_Plane_TS_l_irk_n_erk.hpp"



bool SWE_Plane_TS_l_irk_n_erk::registerShacks(
		sweet::ShackDictionary *io_shackDict
)
{
	ts_l_irk.registerShacks(io_shackDict);
	PDESWEPlaneTS_BaseInterface::registerShacks(io_shackDict);

	return true;
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_irk_n_erk::euler_timestep_update_nonlinear(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t	///< time updates
)
{
	/*
	 * non-conservative (advective) formulation:
	 *
	 *	h_t = -(u*h)_x - (v*h)_y
	 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
	 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
	 */
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);


}


void SWE_Plane_TS_l_irk_n_erk::run_timestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_irk_n_erk: Only constant time step size allowed");

	sweet::PlaneData_Spectral h_linear_t1 = io_h;
	sweet::PlaneData_Spectral u_linear_t1 = io_u;
	sweet::PlaneData_Spectral v_linear_t1 = io_v;

	ts_l_irk.run_timestep(
			h_linear_t1, u_linear_t1, v_linear_t1,
			i_dt,
			i_simulation_timestamp
		);

	// compute non-linear tendencies at half time step
	sweet::PlaneData_Spectral h_dt_nonlinear(ops->planeDataConfig);
	sweet::PlaneData_Spectral u_dt_nonlinear(ops->planeDataConfig);
	sweet::PlaneData_Spectral v_dt_nonlinear(ops->planeDataConfig);

	// standard time stepping
	euler_timestep_update_nonlinear(
			io_h, io_u, io_v,
			h_dt_nonlinear, u_dt_nonlinear, v_dt_nonlinear
		);

	io_h = h_linear_t1 + h_dt_nonlinear*i_dt;
	io_u = u_linear_t1 + u_dt_nonlinear*i_dt;
	io_v = v_linear_t1 + v_dt_nonlinear*i_dt;
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_irk_n_erk::setup(
	sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	assert(io_ops != nullptr);
	assert(shackPDESWETimeDisc != nullptr);
	assert(shackPDESWEPlane != nullptr);

	timestepping_order_linear = shackPDESWETimeDisc->timestepping_order;
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_irk.setup(io_ops, timestepping_order_linear);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_irk_n_erk");


	if (timestepping_order_linear != 1)
		SWEETError("SWE_Plane_TS_l_irk_n_erk: Only 1st order TS supported with this implementation. Please set --timestepping-order 1.");

	timestepping_order_nonlinear = timestepping_order_linear;
	timestepping_rk.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear);
	return true;
}

