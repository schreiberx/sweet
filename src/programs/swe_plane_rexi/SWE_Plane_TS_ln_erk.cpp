/*
 * SWE_Plane_TS_ln_erk.cpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane_rexi.cpp
 *					which was also written by Pedro Peixoto
 */

#include "SWE_Plane_TS_ln_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void SWE_Plane_TS_ln_erk::euler_timestep_update(
		const PlaneData &i_h,	///< prognostic variables (perturbed part of height)
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_h_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double &o_dt,			///< time step restriction
		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	/*
	 * TIME STEP SIZE
	 */
	if (i_fixed_dt <= 0)
		FatalError("Only fixed time step size allowed");

	o_dt = i_fixed_dt;

	// A-grid method
	if (!simVars.disc.use_staggering)
	{
		/*
		 * non-conservative (advective) formulation:
		 *
		 *	h_t = -(u*h)_x - (v*h)_y
		 *	u_t = -g * h_x - u * u_x - v * u_y + f*v
		 *	v_t = -g * h_y - u * v_x - v * v_y - f*u
		 */
		PlaneData total_h = i_h + simVars.sim.h0;

		o_u_t = -simVars.sim.gravitation*op.diff_c_x(total_h) - i_u*op.diff_c_x(i_u) - i_v*op.diff_c_y(i_u);
		o_v_t = -simVars.sim.gravitation*op.diff_c_y(total_h) - i_u*op.diff_c_x(i_v) - i_v*op.diff_c_y(i_v);

		o_u_t += simVars.sim.f0*i_v;
		o_v_t -= simVars.sim.f0*i_u;

		// standard update
		o_h_t = -op.diff_c_x(i_u*total_h) - op.diff_c_y(i_v*total_h);
	}
	else // simVars.disc.use_staggering = true
	{
		FatalError("TODO (MaS@Pedro): The A-grid formulation above assumes that h is the perturbation.");
		// STAGGERED GRID

		PlaneData U(i_h.planeDataConfig);
		PlaneData V(i_h.planeDataConfig);
		PlaneData H(i_h.planeDataConfig);

		/*
		 * Sadourny energy conserving scheme
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
		if (simVars.pde.use_nonlinear_equations > 0) // nonlinear case
		{
			U = op.avg_b_x(i_h)*i_u;
			V = op.avg_b_y(i_h)*i_v;
		}
		else // linear case
		{
			U = simVars.sim.h0*i_u;
			V = simVars.sim.h0*i_v;
		}

		if (simVars.pde.use_nonlinear_equations > 0) // nonlinear case
			H = simVars.sim.gravitation*i_h + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));
		else // linear case
			H = simVars.sim.gravitation*i_h;// + 0.5*(op.avg_f_x(i_u*i_u) + op.avg_f_y(i_v*i_v));


		if (simVars.pde.use_nonlinear_equations > 0) //nonlinear case
		{
			// Potential vorticity
			if (op.avg_b_x(op.avg_b_y(i_h)).reduce_min() < 0.00000001)
			{
				std::cerr << "Test case not adequate for vector invariant formulation. Null or negative water height" << std::endl;
				exit(1);
			}

			PlaneData q = (op.diff_b_x(i_v) - op.diff_b_y(i_u) + simVars.sim.f0) / op.avg_b_x(op.avg_b_y(i_h));

			// u, v tendencies
			// Energy conserving scheme
			o_u_t = op.avg_f_y(q*op.avg_b_x(V)) - op.diff_b_x(H);
			o_v_t = -op.avg_f_x(q*op.avg_b_y(U)) - op.diff_b_y(H);
		}
		else // linear case
		{
			o_u_t = op.avg_f_y(simVars.sim.f0*op.avg_b_x(i_v)) - op.diff_b_x(H);
			o_v_t = -op.avg_f_x(simVars.sim.f0*op.avg_b_y(i_u)) - op.diff_b_y(H);
		}

		/*
		 * P UPDATE
		 */
		if (simVars.pde.use_nonlinear_equations > 0){ //full nonlinear divergence
			// standard update
			o_h_t = -op.diff_f_x(U) - op.diff_f_y(V);
		}
		else // use linear divergence
		{
			o_h_t = -op.diff_f_x(simVars.sim.h0*i_u) - op.diff_f_y(simVars.sim.h0*i_v);
		}
	}
}



void SWE_Plane_TS_ln_erk::run_timestep(
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
		FatalError("SWE_Plane_TS_ln_erk: Only constant time step size allowed");

	if (i_simulation_timestamp + i_fixed_dt > i_max_simulation_time)
		i_fixed_dt = i_max_simulation_time-i_simulation_timestamp;

	o_dt = i_fixed_dt;

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&SWE_Plane_TS_ln_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_h, io_u, io_v,
			o_dt,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp,
			i_max_simulation_time
		);
}



/*
 * Setup
 */
void SWE_Plane_TS_ln_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


SWE_Plane_TS_ln_erk::SWE_Plane_TS_ln_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



SWE_Plane_TS_ln_erk::~SWE_Plane_TS_ln_erk()
{
}

