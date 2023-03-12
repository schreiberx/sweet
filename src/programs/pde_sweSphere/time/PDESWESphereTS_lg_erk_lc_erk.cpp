/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_erk_lc_erk.hpp"



bool PDESWESphereTS_lg_erk_lc_erk::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(io_ops, shackPDESWETimeDisc->timestepping_order);
}


bool PDESWESphereTS_lg_erk_lc_erk::setup_main(
		sweet::SphereOperators *io_ops,
		int i_order	///< order of RK time stepping method
)
{
	ops = io_ops;

	timestepping_order = i_order;
	setupFG();

	return true;
}


void PDESWESphereTS_lg_erk_lc_erk::runTimestep(
		sweet::SphereData_Spectral &io_phi_pert,	///< prognostic variables
		sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.runTimestep(
				this,
				&PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lc,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.runTimestep(
				this,
				&PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lg,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else
	{
		SWEETError("This time stepping order is not yet supported!");
	}
}





/*
 * Main routine for method to be used in case of finite differences
 */
void PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update(
		const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
		const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
		const sweet::SphereData_Spectral &i_div,	///< prognostic variables

		sweet::SphereData_Spectral &o_phi_t,	///< time updates
		sweet::SphereData_Spectral &o_vort_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	if (shackPDESWESphere->sphere_use_fsphere)
	{
		double gh0 = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		o_phi_t = -gh0*i_div;
		o_div_t = -ops->laplace(i_phi);

		o_vort_t = -shackPDESWESphere->sphere_fsphere_f0*i_div;
		o_div_t += shackPDESWESphere->sphere_fsphere_f0*i_vort;
	}
	else
	{
#if 0
		double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		o_phi_t = -gh*i_div;
		o_div_t = -ops->laplace(i_phi);

		/*
		 * This doesn't converge to the reference solution
		 */
		sweet::SphereData_Spectral f(fg);
		o_vort_t = -f*i_div;
		o_div_t += f*i_vort;

#else

		double gh = shackPDESWESphere->gravitation * shackPDESWESphere->h0;

		/*
		 * Apply Coriolis Effect in physical VELOCITY space
		 */
		sweet::SphereData_Physical ug(i_phi.sphereDataConfig);
		sweet::SphereData_Physical vg(i_phi.sphereDataConfig);
		ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

		sweet::SphereData_Physical tmpg1 = ug*fg;
		sweet::SphereData_Physical tmpg2 = vg*fg;

		ops->uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

		o_vort_t *= -1.0;
		o_div_t += -ops->laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh)*i_div;
#endif
	}
}




void PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lg(
		const sweet::SphereData_Spectral &i_phi,	///< prognostic variables
		const sweet::SphereData_Spectral &i_vort,	///< prognostic variables
		const sweet::SphereData_Spectral &i_div,	///< prognostic variables

		sweet::SphereData_Spectral &o_phi_t,	///< time updates
		sweet::SphereData_Spectral &o_vort_t,	///< time updates
		sweet::SphereData_Spectral &o_div_t,	///< time updates

		double i_simulation_timestamp
)
{
	double gh0 = shackPDESWESphere->gravitation*shackPDESWESphere->h0;


	/*
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 * Section "lg_erk"
	 */


	/*
	 * step 1a
	 */
	o_vort_t.spectral_set_zero();

	/*
	 * step 1b
	 */
	o_div_t = -ops->laplace(i_phi);

	/*
	 * step 2a
	 */
	o_phi_t = -gh0*i_div;
}



void PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lc(
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
	 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
	 * Section "lc_erk"
	 */

	/*
	 * step 1a
	 */
	sweet::SphereData_Physical ug(i_phi.sphereDataConfig);
	sweet::SphereData_Physical vg(i_phi.sphereDataConfig);
	ops->vrtdiv_to_uv(i_vort, i_div, ug, vg);

	/*
	 * step 1b
	 */
	sweet::SphereData_Physical tmp_u = ug*fg;
	sweet::SphereData_Physical tmp_v = vg*fg;

	/*
	 * step 1c
	 */
	ops->uv_to_vrtdiv(tmp_u, tmp_v, o_div_t, o_vort_t);

	/*
	 * step 1d
	 */
	o_vort_t *= -1.0;


	/*
	 * step 1e
	 * Nothing to do
	 */

	/*
	 * step 2a
	 * Zero tendencies
	 */
	o_phi_t.spectral_set_zero();

}



/**
 * This routine is used by other time step implementations
 */
void PDESWESphereTS_lg_erk_lc_erk::euler_timestep_update_lc(
		sweet::SphereData_Spectral &io_phi,		///< prognostic variables
		sweet::SphereData_Spectral &io_vort,	///< prognostic variables
		sweet::SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	sweet::SphereData_Spectral tmp_phi(io_phi.sphereDataConfig);
	sweet::SphereData_Spectral tmp_vort(io_vort.sphereDataConfig);
	sweet::SphereData_Spectral tmp_div(io_div.sphereDataConfig);

	euler_timestep_update_lc(
			io_phi,
			io_vort,
			io_div,

			tmp_phi,
			tmp_vort,
			tmp_div,

			i_simulation_timestamp
		);

	io_phi += i_dt*tmp_phi;
	io_vort += i_dt*tmp_vort;
	io_div += i_dt*tmp_div;
}



PDESWESphereTS_lg_erk_lc_erk::PDESWESphereTS_lg_erk_lc_erk()
{
}



PDESWESphereTS_lg_erk_lc_erk::~PDESWESphereTS_lg_erk_lc_erk()
{
}

