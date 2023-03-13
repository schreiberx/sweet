/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_ln_erk.hpp"


bool PDESWESphereTS_ln_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphereTS_ln_erk::setup(
		sweet::SphereOperators *io_ops,
		int i_order	///< order of RK time stepping method
)
{
	ops = io_ops,
	timestepping_order = i_order;

	setupFG();

	return true;
}


/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphereTS_ln_erk::euler_timestep_update_pert(
		const sweet::SphereData_Spectral &i_phi_pert,	///< prognostic variables
		const sweet::SphereData_Spectral &i_vrt,	///< prognostic variables
		const sweet::SphereData_Spectral &i_div,	///< prognostic variables

		sweet::SphereData_Spectral &o_phi_pert_t,	///< time updates
		sweet::SphereData_Spectral &o_vrt_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

	/*
	 * NON-LINEAR SWE
	 *
	 * See
	 * 	Williamson, David L., Drake, John B., Hack, James J., Jakob, Rudiger, & Swarztrauber, Paul N. (1992).
	 * 	A standard test set for numerical approximations to the shallow water equations in spherical geometry.
	 * 	Journal of Computational Physics, 102(1), 211–224. https://doi.org/10.1016/S0021-9991(05)80016-6
	 *
	 * "2.3 Vorticity/Divergence Form"
	 */

	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 */
	sweet::SphereData_Physical phi_pert_phys = i_phi_pert.toPhys();

	/*
	 * Step 1a
	 */
	sweet::SphereData_Physical ug, vg;
	ops->vrtdiv_to_uv(i_vrt, i_div, ug, vg);

	/*
	 * Step 1b
	 */
	sweet::SphereData_Physical vrtg = i_vrt.toPhys();

	/*
	 * Step 1c
	 */

	using namespace sweet;

	// left part of eq. (19)
	sweet::SphereData_Physical u_nl = ug*(vrtg+fg);

	// left part of eq. (20)
	sweet::SphereData_Physical v_nl = vg*(vrtg+fg);

	/*
	 * Step 1d
	 */
	// Eq. (21) & left part of Eq. (22)
	ops->uv_to_vrtdiv(u_nl, v_nl, o_div_t, o_vrt_t);


	/*
	 * Step 1e
	 */
	o_vrt_t *= -1.0;

	/*
	 * Step 1f
	 */
	// Right part of Eq. (22)
	sweet::SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	sweet::SphereData_Spectral e = phi_pert_phys+tmpg;

	/*
	 * Step 1g
	 */
	o_div_t -= ops->laplace(e);

	/*
	 * Compute Phi geopotential tendencies
	 */

	/*
	 * Step 2a
	 */
	u_nl = ug*(phi_pert_phys + gh0);
	v_nl = vg*(phi_pert_phys + gh0);

	ops->uv_to_vrtdiv(u_nl,v_nl, e, o_phi_pert_t);

	o_phi_pert_t *= -1.0;


}



void PDESWESphereTS_ln_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi,		///< prognostic variables
		sweet::SphereData_Spectral &io_vort,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,		///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		SWEETError("Only constant time step size allowed");

	// standard time stepping
	timestepping_rk.runTimestep(
			this,
			&PDESWESphereTS_ln_erk::euler_timestep_update_pert,	///< pointer to function to compute euler time step updates
			io_phi, io_vort, io_div,
			i_fixed_dt,
			timestepping_order,
			i_simulation_timestamp
		);
}


PDESWESphereTS_ln_erk::PDESWESphereTS_ln_erk()
{
}



PDESWESphereTS_ln_erk::~PDESWESphereTS_ln_erk()
{
}

