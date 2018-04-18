/*
 * Adv_Plane_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include "Adv_Plane_TS_na_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void Adv_Plane_TS_na_erk::euler_timestep_update(
		const PlaneData &i_phi,	///< prognostic variables
		const PlaneData &i_u,	///< prognostic variables
		const PlaneData &i_v,	///< prognostic variables

		PlaneData &o_phi_t,	///< time updates
		PlaneData &o_u_t,	///< time updates
		PlaneData &o_v_t,	///< time updates

		double i_simulation_timestamp
)
{
	/**
	 * We simply compute
	 * 	-DIV(rho*U) = -rho DIV(U) - U.GRAD(rho) = - U.GRAD(rho)
	 * which is the Lagrangian contribution only.
	 *
	 * This is the case because the velocity field is divergence free!!!
	 */
/*
	PlaneData* ext_u;
	PlaneData* ext_v;

	simVars.sim.getExternalForces(0, &ext_u);
	simVars.sim.getExternalForces(1, &ext_v);
*/
	o_phi_t = -op.diff_c_x(i_phi*i_u) - op.diff_c_y(i_phi*i_v);

	o_u_t.spectral_set_zero();
	o_v_t.spectral_set_zero();
}



void Adv_Plane_TS_na_erk::run_timestep(
		PlaneData &io_phi,		///< prognostic variables
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&Adv_Plane_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_u, io_v,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}



/*
 * Setup
 */
void Adv_Plane_TS_na_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


Adv_Plane_TS_na_erk::Adv_Plane_TS_na_erk(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Adv_Plane_TS_na_erk::~Adv_Plane_TS_na_erk()
{
}

