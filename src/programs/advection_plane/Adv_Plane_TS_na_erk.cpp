/*
 * Adv_Plane_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "Adv_Plane_TS_na_erk.hpp"





/*
 * Main routine for method to be used in case of finite differences
 */
void Adv_Plane_TS_na_erk::euler_timestep_update(
		const PlaneData_Spectral &i_phi,	///< prognostic variables
		const PlaneData_Spectral &i_u,	///< prognostic variables
		const PlaneData_Spectral &i_v,	///< prognostic variables

		PlaneData_Spectral &o_phi_t,	///< time updates
		PlaneData_Spectral &o_u_t,	///< time updates
		PlaneData_Spectral &o_v_t,	///< time updates

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

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		PlaneData_Spectral u(i_phi.planeDataConfig);
		PlaneData_Spectral v(i_phi.planeDataConfig);

		simVars.benchmark.getExternalForcesCallback(1, simVars.timecontrol.current_simulation_time, &u, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, simVars.timecontrol.current_simulation_time, &v, simVars.benchmark.getExternalForcesUserData);

		o_phi_t = -op.diff_c_x(i_phi*u) - op.diff_c_y(i_phi*v);
	}
	else
	{
		o_phi_t = -op.diff_c_x(i_phi*i_u) - op.diff_c_y(i_phi*i_v);
	}

//#if SWEET_USE_PLANE_SPECTRAL_SPACE
	o_u_t.spectral_set_zero();
	o_v_t.spectral_set_zero();
///#else
///	o_u_t.physical_set_zero();
///	o_v_t.physical_set_zero();
///#endif
}



void Adv_Plane_TS_na_erk::run_timestep(
		PlaneData_Spectral &io_phi,		///< prognostic variables
		PlaneData_Spectral &io_u,	///< prognostic variables
		PlaneData_Spectral &io_v,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	if (i_dt <= 0)
		SWEETError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&Adv_Plane_TS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (simVars.benchmark.getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		simVars.benchmark.getExternalForcesCallback(1, simVars.timecontrol.current_simulation_time+i_dt, &io_u, simVars.benchmark.getExternalForcesUserData);
		simVars.benchmark.getExternalForcesCallback(2, simVars.timecontrol.current_simulation_time+i_dt, &io_v, simVars.benchmark.getExternalForcesUserData);
	}
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

