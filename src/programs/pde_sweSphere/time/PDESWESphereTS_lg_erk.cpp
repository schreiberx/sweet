/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_erk.hpp"




bool PDESWESphereTS_lg_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphereTS_lg_erk::setup_main(
		const sweet::SphereOperators *io_ops,
		int i_order	///< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_order;
	return true;
}




void PDESWESphereTS_lg_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	// standard time stepping
	timestepping_rk.runTimestep(
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
		const sweet::SphereData_Spectral &i_phi,
		const sweet::SphereData_Spectral &i_vort,
		const sweet::SphereData_Spectral &i_div,

		sweet::SphereData_Spectral &o_phi_t,	///< time updates
		sweet::SphereData_Spectral &o_vort_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	o_phi_t = -gh*i_div;
	o_div_t = -ops->laplace(i_phi);
	o_vort_t.spectral_set_zero();
}



PDESWESphereTS_lg_erk::PDESWESphereTS_lg_erk()
{
}


PDESWESphereTS_lg_erk::~PDESWESphereTS_lg_erk()
{
}

