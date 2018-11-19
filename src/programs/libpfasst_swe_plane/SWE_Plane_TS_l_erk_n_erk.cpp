/*
 * SWE_Plane_TS_l_erk_n_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "../swe_plane/SWE_Plane_TS_l_erk_n_erk.hpp"

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
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * linearized non-conservative (advective) formulation:
	 *
	 * h_t = -h0*u_x - h0*v_ym
	 * u_t = -g * h_x + f*v
	 * v_t = -g * h_y - f*u
	 */
	o_h_t = -(op.diff_c_x(i_u) + op.diff_c_y(i_v))*simVars.sim.h0;
	o_u_t = -simVars.sim.gravitation*op.diff_c_x(i_h) + simVars.sim.plane_rotating_f0*i_v;
	o_v_t = -simVars.sim.gravitation*op.diff_c_y(i_h) - simVars.sim.plane_rotating_f0*i_u;

}



void SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

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
	//o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
	o_u_t = -i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
	o_v_t = -i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);
	if (simVars.pde.use_only_linear_divergence == 1) //only nonlinear advection left to solve
		o_h_t = - (i_u*op.diff_c_x(i_h) + i_v*op.diff_c_y(i_h));
	else //full nonlinear equation on h
		o_h_t = -op.diff_c_x(i_u*i_h) - op.diff_c_y(i_v*i_h);
}



void SWE_Plane_TS_l_erk_n_erk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
			FatalError("Only fixed time step size allowed (set --dt)");

	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				1,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
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
				io_h, io_u, io_v,
				i_dt*0.5,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Plane_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_h, io_u, io_v,
				i_dt*0.5,
				2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		FatalError("SWE_Plane_TS_l_erk_n_erk: Order not yet supported!");
	}
}





/*
 * Setup
 */
void SWE_Plane_TS_l_erk_n_erk::setup(
		int i_order,	///< order of RK time stepping method
		int i_order2	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;

	if (timestepping_order != timestepping_order2)
		FatalError("TODO: Currently, both time stepping orders (1 and 2) have to match!");

	timestepping_rk_linear.setupBuffers(op.planeDataConfig, timestepping_order);
	timestepping_rk_nonlinear.setupBuffers(op.planeDataConfig, timestepping_order2);

	if (simVars.disc.use_staggering)
		FatalError("SWE_Plane_TS_l_erk_n_erk: Staggering not supported");
}



SWE_Plane_TS_l_erk_n_erk::SWE_Plane_TS_l_erk_n_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order, simVars.disc.timestepping_order2);
}



SWE_Plane_TS_l_erk_n_erk::~SWE_Plane_TS_l_erk_n_erk()
{
}

