/*
 * SWE_Plane_TS_l_rexi_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_rexi_n_erk.hpp"

#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


bool SWE_Plane_TS_l_rexi_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	ts_l_rexi.shackRegistration(io_shackDict);
	PDESWEPlaneTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_rexi_n_erk::euler_timestep_update_nonlinear(
		const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

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
	//o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);
	o_u_t = -i_u*ops->diff_c_x(i_u) - i_v*ops->diff_c_y(i_u);
	o_v_t = -i_u*ops->diff_c_x(i_v) - i_v*ops->diff_c_y(i_v);
	if (use_only_linear_divergence) //only nonlinear advection left to solve
		o_h_t = - (i_u*ops->diff_c_x(i_h) + i_v*ops->diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -ops->diff_c_x(i_u*i_h) - ops->diff_c_y(i_v*i_h);

}


void SWE_Plane_TS_l_rexi_n_erk::run_timestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Only constant time step size allowed");

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
		SWEETError("SWE_Plane_TS_l_rexi_n_erk: Explicit erk order not implemented for this scheme, please set --timestepping-order2 to 1 or 2.");
	}
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_rexi_n_erk::setup(
	sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	ts_l_rexi.setup(io_ops, "phi0");

	timestepping_order_nonlinear = shackPDESWETimeDisc->timestepping_order;
	timestepping_rk.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("Staggering not supported for l_rexi_n_erk");

	return true;
}

