/*
 * Adv_Plane_TS_ln_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <programs/pdeAdvectionPlane/PDEAdvPlaneTS_na_erk.hpp>



/*
 * Main routine for method to be used in case of finite differences
 */
void PDEAdvPlaneTS_na_erk::euler_timestep_update(
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

	if (shackBenchmark->getExternalForcesCallback != nullptr)
	{
		PlaneData_Spectral u(i_phi.planeDataConfig);
		PlaneData_Spectral v(i_phi.planeDataConfig);

		shackBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->current_simulation_time,
				&u,
				shackBenchmark
			);
		shackBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->current_simulation_time,
				&v,
				shackBenchmark
			);

		o_phi_t = -op.diff_c_x(i_phi*u) - op.diff_c_y(i_phi*v);
	}
	else
	{
		o_phi_t = -op.diff_c_x(i_phi*i_u) - op.diff_c_y(i_phi*i_v);
	}

	o_u_t.spectral_set_zero();
	o_v_t.spectral_set_zero();
}


void PDEAdvPlaneTS_na_erk::run_timestep(
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
			&PDEAdvPlaneTS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (shackBenchmark->getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		shackBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->current_simulation_time+i_dt,
				&io_u,
				shackBenchmark
			);
		shackBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->current_simulation_time+i_dt,
				&io_v,
				shackBenchmark
			);
	}
}



/*
 * Setup
 */
void PDEAdvPlaneTS_na_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


PDEAdvPlaneTS_na_erk::PDEAdvPlaneTS_na_erk(
		sweet::ShackDictionary &io_shackDict,
		PlaneOperators &i_op
)	:
		op(i_op)
{
	shackTimeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();
	shackTimestepControl = io_shackDict.getAutoRegistration<ShackTimestepControl>();
	shackBenchmark = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();

	setup(shackTimeDisc->timestepping_order);
}



PDEAdvPlaneTS_na_erk::~PDEAdvPlaneTS_na_erk()
{
}

