/*
 * PDESWESphereTS_lg_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_erk.hpp"



void PDESWESphereTS_lg_erk::run_timestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.run_timestep(
			this,
			&PDESWESphereTS_lg_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}

/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphereTS_lg_erk::euler_timestep_update(
		const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
		const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
		const sweet::SphereData_Spectral &i_div,	///< prognostic variables

		sweet::SphereData_Spectral &o_phi_t,	///< time updates
		sweet::SphereData_Spectral &o_vort_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = shackDict.sim.gravitation * shackDict.sim.h0;

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vort_t.spectral_set_zero();
}





/*
 * Setup
 */
void PDESWESphereTS_lg_erk::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;
}


PDESWESphereTS_lg_erk::PDESWESphereTS_lg_erk(
		sweet::ShackDictionary &i_shackDict,
		sweet::SphereOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op)
{
	setup(timestepping_order);
}



PDESWESphereTS_lg_erk::~PDESWESphereTS_lg_erk()
{
}

