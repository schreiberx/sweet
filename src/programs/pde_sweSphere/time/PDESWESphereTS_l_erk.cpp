/*
 * PDESWESphereTS_l_erk.cpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_l_erk.hpp"

bool PDESWESphereTS_l_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphereTS_l_erk::setup_main(
		sweet::SphereOperators *io_ops,
		int i_order	///< order of RK time stepping method
)
{
	ops = io_ops;

	setupFG();

	timestepping_order = i_order;
	return true;
}



void PDESWESphereTS_l_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphereTS_l_erk::euler_timestep_update,	///< pointer to function to compute euler time step updates
			io_phi_pert, io_vrt, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphereTS_l_erk::euler_timestep_update(
		const sweet::SphereData_Spectral &i_phi_pert,	///< prognostic variables
		const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
		const sweet::SphereData_Spectral &i_div,	///< prognostic variables

		sweet::SphereData_Spectral &o_phi_pert_t,	///< time updates
		sweet::SphereData_Spectral &o_vort_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (!shackPDESWESphere->sphere_use_fsphere)
	{

		double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
		 */
		/*
		 * Step 1a
		 */
		sweet::SphereData_Physical ug(i_phi_pert.sphereDataConfig);
		sweet::SphereData_Physical vg(i_phi_pert.sphereDataConfig);
		ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		/*
		 * Step 1b
		 */
		sweet::SphereData_Physical tmpg1 = ug*fg;
		sweet::SphereData_Physical tmpg2 = vg*fg;

		/*
		 * Step 1c
		 */
		ops->uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		/*
		 * Step 1d
		 */
		o_vort_t *= -1.0;

		/*
		 * Step 1e
		 */
		o_div_t += -ops->laplace(i_phi_pert);

		/*
		 * DIV on velocity field
		 */
		o_phi_pert_t = (-gh0)*i_div;
	}
	else
	{
		double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		o_div_t = -ops->laplace(i_phi_pert);

		o_vort_t = -shackPDESWESphere->sphere_fsphere_f0*i_div;
		o_div_t += shackPDESWESphere->sphere_fsphere_f0*i_vort;

		o_phi_pert_t = -gh*i_div;
	}
}




PDESWESphereTS_l_erk::PDESWESphereTS_l_erk()
{
}



PDESWESphereTS_l_erk::~PDESWESphereTS_l_erk()
{
}

