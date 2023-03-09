/*
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDEAdvectionPlaneTS_na_erk.hpp"


PDEAdvectionPlaneTS_na_erk::PDEAdvectionPlaneTS_na_erk()
{
}


PDEAdvectionPlaneTS_na_erk::~PDEAdvectionPlaneTS_na_erk()
{
}


/*
 * Setup
 */
bool PDEAdvectionPlaneTS_na_erk::setup(sweet::PlaneOperators *io_ops)
{
	PDEAdvectionPlaneTS_BaseInterface::setup(io_ops);
	timestepping_order = shackPDEAdvTimeDisc->timestepping_order;
	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDEAdvectionPlaneTS_na_erk::euler_timestep_update(
		const sweet::PlaneData_Spectral &i_phi,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
		const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

		sweet::PlaneData_Spectral &o_phi_t,	///< time updates
		sweet::PlaneData_Spectral &o_u_t,	///< time updates
		sweet::PlaneData_Spectral &o_v_t,	///< time updates

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

	if (shackPDEAdvBenchmark->getExternalForcesCallback != nullptr)
	{
		sweet::PlaneData_Spectral u(i_phi.planeDataConfig);
		sweet::PlaneData_Spectral v(i_phi.planeDataConfig);

		shackPDEAdvBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->current_simulation_time,
				&u,
				shackPDEAdvBenchmark
			);
		shackPDEAdvBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->current_simulation_time,
				&v,
				shackPDEAdvBenchmark
			);

		o_phi_t = -ops->diff_c_x(i_phi*u) - ops->diff_c_y(i_phi*v);
	}
	else
	{
		o_phi_t = -ops->diff_c_x(i_phi*i_u) - ops->diff_c_y(i_phi*i_v);
	}

	o_u_t.spectral_set_zero();
	o_v_t.spectral_set_zero();
}


void PDEAdvectionPlaneTS_na_erk::runTimestep(
		sweet::PlaneData_Spectral &io_phi,		///< prognostic variables
		sweet::PlaneData_Spectral &io_u,	///< prognostic variables
		sweet::PlaneData_Spectral &io_v,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	assert(i_dt > 0);

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDEAdvectionPlaneTS_na_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi, io_u, io_v,
			i_dt,
			timestepping_order,
			i_simulation_timestamp
		);

	if (shackPDEAdvBenchmark->getExternalForcesCallback != nullptr)
	{
		// this is just called for cosmetic reasons to update the velocity field
		shackPDEAdvBenchmark->getExternalForcesCallback(
				1,
				shackTimestepControl->current_simulation_time+i_dt,
				&io_u,
				shackPDEAdvBenchmark
			);
		shackPDEAdvBenchmark->getExternalForcesCallback(
				2,
				shackTimestepControl->current_simulation_time+i_dt,
				&io_v,
				shackPDEAdvBenchmark
			);
	}
}


