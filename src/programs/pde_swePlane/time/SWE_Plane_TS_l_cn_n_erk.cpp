/*
 * SWE_Plane_TS_l_cn_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_cn_n_erk.hpp"

#include <sweet/core/plane/PlaneOperatorsComplex.hpp>


bool SWE_Plane_TS_l_cn_n_erk::shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{

	ts_l_cn.shackRegistration(io_shackDict);
	PDESWEPlaneTS_BaseInterface::shackRegistration(io_shackDict);

	return true;
}

/*
 *
 * U_t=LU+N(U)
 *
 * Strang Splitting used

 *
 * 2nd order scheme
 * 1) dt/2 step of linear part with CN
 * 2) dt step of nonlinear part with RK from step 1)
 * 3) dt/2 step of linear part with CN from step 2)
 * Formally: U(n+1)=exp(dtL/2)exp(dtN)exp(dtL/2)U(n)
 *
 */

void SWE_Plane_TS_l_cn_n_erk::euler_timestep_update_nonlinear(
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


void SWE_Plane_TS_l_cn_n_erk::run_timestep(
		sweet::PlaneData_Spectral &io_h,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("SWE_Plane_TS_l_cn_n_erk: Only constant time step size allowed (set --dt)");

	// Half one for linear parts
	ts_l_cn.run_timestep(
			io_h, io_u, io_v,
			i_dt*0.5,
			i_simulation_timestamp
		);

	// standard time stepping for nonlinear parts
	timestepping_rk.run_timestep(
			this,
			&SWE_Plane_TS_l_cn_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order_nonlinear,
			i_simulation_timestamp
		);

	// Half one for linear parts
	ts_l_cn.run_timestep(
			io_h, io_u, io_v,
			i_dt*0.5,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
bool SWE_Plane_TS_l_cn_n_erk::setup(
		sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);

	timestepping_order_linear = shackPDESWETimeDisc->timestepping_order;
	timestepping_order_nonlinear = shackPDESWETimeDisc->timestepping_order2;
	crank_nicolson_damping_factor = 0.5;
	use_only_linear_divergence = shackPDESWEPlane->use_only_linear_divergence;

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_cn_n_erk: Staggering not supported for l_cn_n_erk");

	if (timestepping_order_linear > 0 && timestepping_order_linear != 2)
		std::cout << "SWE_Plane_TS_l_cn_n_erk Warning: Ignoring timestepping_order for linear time integration because this is always given by the Crank-Nicolson scheme" <<std::endl;

	if (timestepping_order_nonlinear <= 0)
		SWEETError("SWE_Plane_TS_l_cn_n_erk: Please set --timestepping-order2 to define the order of the nonlinear time integration");

	//ts_l_cn.setup(2, i_crank_nicolson_damping_factor);
	ts_l_cn.setup(io_ops);

	timestepping_rk.setupBuffers(ops->planeDataConfig, timestepping_order_nonlinear);
	return true;
}

