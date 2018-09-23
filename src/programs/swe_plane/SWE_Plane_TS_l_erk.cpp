/*
 * SWE_Plane_TS_l_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_l_erk.hpp"


/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_l_erk::euler_timestep_update(
		const PlaneData &i_h,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	// A- grid method
	if (!simVars.disc.use_staggering)
	{
		/*
		 * linearized non-conservative (advective) formulation:
		 *
		 * h_t = -h0*u_x - h0*v_ym
		 * u_t = -g * h_x + f*v
		 * v_t = -g * h_y - f*u
		 */

#if 1
		o_u_t = -simVars.sim.gravitation*op.diff_c_x(i_h) + simVars.sim.f0*i_v;
		o_v_t = -simVars.sim.gravitation*op.diff_c_y(i_h) - simVars.sim.f0*i_u;

		// standard update
		o_h_t = -(op.diff_c_x(i_u) + op.diff_c_y(i_v))*simVars.sim.h0;
#else

	#if 0
		// U-only
		o_u_t = -simVars.sim.gravitation*op.diff_c_x(i_h) + simVars.sim.f0*i_v;
		//o_v_t.physical_set_zero();
		o_v_t = - simVars.sim.f0*i_u;

		// standard update
		o_h_t = -(op.diff_c_x(i_u))*simVars.sim.h0;

	#else
		// V-only
		//o_u_t.spectral_set_zero();
		o_u_t = +simVars.sim.f0*i_v;
		o_v_t = -simVars.sim.gravitation*op.diff_c_y(i_h) - simVars.sim.f0*i_u;// - simVars.sim.f0*i_u;

		// standard update
		o_h_t = -(op.diff_c_y(i_v))*simVars.sim.h0;
	#endif

#endif
	}
	else // simVars.disc.use_staggering = true
	{
		// STAGGERED GRID

		PlaneData U(i_h.planeDataConfig);
		PlaneData V(i_h.planeDataConfig);
		PlaneData H(i_h.planeDataConfig);

		/* Sadourny energy conserving scheme
		 *
		 * Note, that this grid does not follow the formulation
		 * in the paper of Robert Sadourny, but looks as follows:
		 *
		 *              ^
		 *              |
		 *       ______v0,1_____
		 *       |             |
		 *       |			   |
		 *       |             |
		 *  u0,0 |->  H/P0,0   |u1,0 ->
		 *(0,0.5)|			   |
		 *       |      ^      |
		 *   q0,0|______|______|
		 * (0,0)      v0,0
		 *           (0.5,0)
		 *
		 * V_t + q N x (P V) + grad( g P + 1/2 V*V) = 0
		 * P_t + div(P V) = 0
		 */

		/*
		 * U and V updates
		 */
		U = simVars.sim.h0*i_u;
		V = simVars.sim.h0*i_v;

		H = simVars.sim.gravitation*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));

		o_u_t = op.avg_f_y(simVars.sim.f0*op.avg_b_x(i_v)) - op.diff_b_x(H);
		o_v_t = -op.avg_f_x(simVars.sim.f0*op.avg_b_y(i_u)) - op.diff_b_y(H);

		/*
		 * P UPDATE
		 */
		o_h_t = -op.diff_f_x(simVars.sim.h0*i_u) - op.diff_f_y(simVars.sim.h0*i_v);
	}
}



void SWE_Plane_TS_l_erk::run_timestep(
		PlaneData &io_h,	///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables

		double i_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		FatalError("SWE_Plane_TS_l_erk: Only constant time step size allowed (please set --dt)");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Plane_TS_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void SWE_Plane_TS_l_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
	timestepping_rk.setupBuffers(op.planeDataConfig, timestepping_order);

	//if (simVars.disc.use_staggering)
	//	FatalError("Staggering not supported for l_erk");
}


SWE_Plane_TS_l_erk::SWE_Plane_TS_l_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Plane_TS_l_erk::~SWE_Plane_TS_l_erk()
{
}

