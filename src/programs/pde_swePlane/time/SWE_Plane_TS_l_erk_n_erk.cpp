/*
 * SWE_Plane_TS_l_erk_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_erk_n_erk.hpp"

/*
 *
 * U_t=LU+N(U)
 *
 * Strang Splitting used
 * 1st order scheme:
 * 1) dt step of linear part with RK1
 * 2) dt step of nonlinear part from step 1) with RK1
 * Formally: U(n+1)=exp(dtN)exp(dtL)U(n)
 *
 * 2nd order scheme
 * 1) dt/2 step of linear part with RK2
 * 2) dt step of nonlinear part with RK2 from step 1)
 * 3) dt/2 step of linear part with RK2 from step 2)
 * Formally: U(n+1)=exp(dtL/2)exp(dtN)exp(dtL/2)U(n)
 *
 */

void SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear(
		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * linearized non-conservative (advective) formulation:
	 *
	 * h_pert_t = -h0*u_x - h0*v_ym
	 * u_t = -g * h_x + f*v
	 * v_t = -g * h_y - f*u
	 */
	o_h_pert_t = -(ops->diff_c_x(i_u) + ops->diff_c_y(i_v))*shackPDESWEPlane->h0;
	o_u_t = -shackPDESWEPlane->gravitation*ops->diff_c_x(i_h_pert) + shackPDESWEPlane->plane_rotating_f0*i_v;
	o_v_t = -shackPDESWEPlane->gravitation*ops->diff_c_y(i_h_pert) - shackPDESWEPlane->plane_rotating_f0*i_u;

}



void SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear(
		const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_h_pert_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

		double i_simulation_timestamp
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
		o_h_pert_t = - (i_u*ops->diff_c_x(i_h_pert) + i_v*ops->diff_c_y(i_h_pert));
	else //full nonlinear equation on h
		o_h_pert_t = -ops->diff_c_x(i_u*i_h_pert) - ops->diff_c_y(i_v*i_h_pert);
}



void SWE_Plane_TS_l_erk_n_erk::run_timestep(
		sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,	///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
			SWEETError("Only fixed time step size allowed (set --dt)");

	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_h_pert, io_u, io_v,
				i_dt,
				1,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert, io_u, io_v,
				i_dt,
				1,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_h_pert, io_u, io_v,
				i_dt*0.5,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h_pert, io_u, io_v,
				i_dt,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_h_pert, io_u, io_v,
				i_dt*0.5,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETError("SWE_Plane_TS_l_erk_n_erk: Order not yet supported!");
	}
}





/*
 * Setup
 */
bool SWE_Plane_TS_l_erk_n_erk::setup(
	sweet::PlaneOperators *io_ops
)
{
	PDESWEPlaneTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;

	if (shackPDESWETimeDisc->timestepping_method == "l_erk_na_ld2_erk")
		use_only_linear_divergence = true;
	else
		use_only_linear_divergence = false;

	if (timestepping_order != timestepping_order2)
		SWEETError("TODO: Currently, both time stepping orders (1 and 2) have to match!");

	timestepping_rk_linear.setupBuffers(ops->planeDataConfig, timestepping_order);
	timestepping_rk_nonlinear.setupBuffers(ops->planeDataConfig, timestepping_order2);

	if (shackPlaneDataOps->space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_erk_n_erk: Staggering not supported");

	return true;
}

