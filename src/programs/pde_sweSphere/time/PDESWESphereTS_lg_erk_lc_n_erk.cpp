/*
 * PDESWESphereTS_lg_erk_lf_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include <sweet/core/StopwatchBox.hpp>
#include "PDESWESphereTS_lg_erk_lc_n_erk.hpp"




void PDESWESphereTS_lg_erk_lc_n_erk::run_timestep(
		sweet::SphereData_Spectral &io_phi_pert,
		sweet::SphereData_Spectral &io_vrt,
		sweet::SphereData_Spectral &io_div,

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{

		if (version_id == 0)
		{
			// HALF time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else if (version_id == 1)
		{
			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);


			// FULL time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);
		}
		else
		{
			SWEETError("Invalid verison id");
		}
	}
	else
	{
		SWEETError("Not yet supported!");
	}
}



void PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_linear(
		const sweet::SphereData_Spectral &i_phi,
		const sweet::SphereData_Spectral &i_vrt,
		const sweet::SphereData_Spectral &i_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vrt_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = shackDict.sim.gravitation * shackDict.sim.h0;

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vrt_t.spectral_set_zero();
}


void PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
		const sweet::SphereData_Spectral &i_phi,
		const sweet::SphereData_Spectral &i_vort,
		const sweet::SphereData_Spectral &i_div,

		sweet::SphereData_Spectral &o_phi_t,
		sweet::SphereData_Spectral &o_vort_t,
		sweet::SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	StopwatchBox::getInstance().main_timestepping_nonlinearities.start();
#endif

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	sweet::SphereData_Physical ug(i_phi.sphereDataConfig);
	sweet::SphereData_Physical vg(i_phi.sphereDataConfig);

	sweet::SphereData_Physical vrtg = i_vort.toPhys();
	sweet::SphereData_Physical divg = i_div.toPhys();
	op.vrtdiv_to_uv(i_vort, i_div, ug, vg);

	sweet::SphereData_Physical phig = i_phi.toPhys();

	sweet::SphereData_Physical tmpg1 = ug*(vrtg+op.fg);
	sweet::SphereData_Physical tmpg2 = vg*(vrtg+op.fg);

	op.uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	sweet::SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	sweet::SphereData_Physical tmpg(i_phi.sphereDataConfig);
	tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = tmpg;
	o_div_t += -op.laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	StopwatchBox::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



/**
 * This routine is used by other time step implementations
 */
void PDESWESphereTS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
		sweet::SphereData_Spectral &io_phi,
		sweet::SphereData_Spectral &io_vort,
		sweet::SphereData_Spectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	sweet::SphereData_Spectral tmp_phi(io_phi.sphereDataConfig);
	sweet::SphereData_Spectral tmp_vort(io_vort.sphereDataConfig);
	sweet::SphereData_Spectral tmp_div(io_div.sphereDataConfig);

	euler_timestep_update_lc_n(
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





/*
 * Setup
 */
void PDESWESphereTS_lg_erk_lc_n_erk::setup(
		int i_order,	///< order of RK time stepping method
		int i_version_id
)
{
	timestepping_order = i_order;

	version_id = i_version_id;
}


PDESWESphereTS_lg_erk_lc_n_erk::PDESWESphereTS_lg_erk_lc_n_erk(
		sweet::ShackDictionary &i_shackDict,
		sweet::SphereOperators &i_op
)	:
		shackDict(i_shackDict),
		op(i_op),
		timestepping_order(-1)
{
}



PDESWESphereTS_lg_erk_lc_n_erk::~PDESWESphereTS_lg_erk_lc_n_erk()
{
}

