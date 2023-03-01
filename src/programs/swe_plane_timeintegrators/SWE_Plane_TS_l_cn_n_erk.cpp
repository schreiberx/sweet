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

#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_cn_n_erk.hpp"

#include <sweet/plane/PlaneOperatorsComplex.hpp>

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
		const PlaneData_Spectral &i_h,	///< prognostic variables
		const PlaneData_Spectral &i_u,	///< prognostic variables
		const PlaneData_Spectral &i_v,	///< prognostic variables

		PlaneData_Spectral &o_h_t,	///< time updates
		PlaneData_Spectral &o_u_t,	///< time updates
		PlaneData_Spectral &o_v_t,	///< time updates

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


void SWE_Plane_TS_l_cn_n_erk::run_timestep(
		PlaneData_Spectral &io_h,	///< prognostic variables
		PlaneData_Spectral &io_u,	///< prognostic variables
		PlaneData_Spectral &io_v,	///< prognostic variables

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
void SWE_Plane_TS_l_cn_n_erk::setup(
		int i_l_order,
		int i_n_order,
		double i_crank_nicolson_damping_factor,
		bool i_use_only_linear_divergence
)
{
	timestepping_order_linear = i_l_order;
	timestepping_order_nonlinear = i_n_order;
	crank_nicolson_damping_factor = i_crank_nicolson_damping_factor;
	use_only_linear_divergence = i_use_only_linear_divergence;

	if (simVars.disc.space_grid_use_c_staggering)
		SWEETError("SWE_Plane_TS_l_cn_n_erk: Staggering not supported for l_cn_n_erk");

	if (timestepping_order_linear > 0 && timestepping_order_linear != 2)
		std::cout<<"SWE_Plane_TS_l_cn_n_erk Warning: Ignoring timestepping_order for linear time integration because this is always given by the Crank-Nicolson scheme"<<std::endl;

	if (timestepping_order_nonlinear <= 0)
		SWEETError("SWE_Plane_TS_l_cn_n_erk: Please set --timestepping-order2 to define the order of the nonlinear time integration");

	//ts_l_cn.setup(2, i_crank_nicolson_damping_factor);
	ts_l_cn.setup(i_crank_nicolson_damping_factor);

	timestepping_rk.setupBuffers(op.planeDataConfig, timestepping_order_nonlinear);
}


SWE_Plane_TS_l_cn_n_erk::SWE_Plane_TS_l_cn_n_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		ts_l_cn(simVars, op)
{
////#if !SWEET_PARAREAL
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2, simVars.disc.timestepping_crank_nicolson_filter, false);
////#endif
}



SWE_Plane_TS_l_cn_n_erk::~SWE_Plane_TS_l_cn_n_erk()
{
}

