/*
 * SWE_Sphere_TS_lg_erk_lf_n_erk.cpp
 *
 *  Created on: 21 Aug 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"

#include <sweet/SimulationBenchmarkTiming.hpp>



void SWE_Sphere_TS_lg_erk_lc_n_erk::run_timestep(
		SphereData_Spectral &io_phi_pert,
		SphereData_Spectral &io_vrt,
		SphereData_Spectral &io_div,

		double i_fixed_dt,			///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
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
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
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
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt*0.5,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);


			// FULL time step for linear part
			timestepping_rk_linear.run_timestep(
					this,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
					io_phi_pert, io_vrt, io_div,
					i_fixed_dt,
					timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
					i_simulation_timestamp
				);

			// HALF time step for non-linear part
			timestepping_rk_nonlinear.run_timestep(
					this,
					&SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n,	///< pointer to function to compute euler time step updates
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



void SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_linear(
		const SphereData_Spectral &i_phi,
		const SphereData_Spectral &i_vrt,
		const SphereData_Spectral &i_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vrt_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
	/*
	 * LINEAR
	 */
	double gh = simVars.sim.gravitation * simVars.sim.h0;

	o_phi_t = -gh*i_div;
	o_div_t = -op.laplace(i_phi);
	o_vrt_t.spectral_set_zero();
}


void SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
		const SphereData_Spectral &i_phi,
		const SphereData_Spectral &i_vort,
		const SphereData_Spectral &i_div,

		SphereData_Spectral &o_phi_t,
		SphereData_Spectral &o_vort_t,
		SphereData_Spectral &o_div_t,

		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.start();
#endif

	/*
	 * NON-LINEAR
	 *
	 * Follows Hack & Jakob formulation
	 */

	SphereData_Physical ug(i_phi.sphereDataConfig);
	SphereData_Physical vg(i_phi.sphereDataConfig);

	SphereData_Physical vrtg = i_vort.toPhys();
	SphereData_Physical divg = i_div.toPhys();
	op.vrtdiv_to_uv(i_vort, i_div, ug, vg);

	SphereData_Physical phig = i_phi.toPhys();

	SphereData_Physical tmpg1 = ug*(vrtg+op.fg);
	SphereData_Physical tmpg2 = vg*(vrtg+op.fg);

	op.uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vort_t);

	o_vort_t *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_t);

	o_phi_t *= -1.0;

	SphereData_Physical tmpg(i_phi.sphereDataConfig);
	tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = tmpg;
	o_div_t += -op.laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



/**
 * This routine is used by other time step implementations
 */
void SWE_Sphere_TS_lg_erk_lc_n_erk::euler_timestep_update_lc_n(
		SphereData_Spectral &io_phi,
		SphereData_Spectral &io_vort,
		SphereData_Spectral &io_div,

		double i_dt,
		double i_simulation_timestamp
)
{
	SphereData_Spectral tmp_phi(io_phi.sphereDataConfig);
	SphereData_Spectral tmp_vort(io_vort.sphereDataConfig);
	SphereData_Spectral tmp_div(io_div.sphereDataConfig);

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
void SWE_Sphere_TS_lg_erk_lc_n_erk::setup(
		int i_order,	///< order of RK time stepping method
		int i_version_id
)
{
	timestepping_order = i_order;

	version_id = i_version_id;
}


SWE_Sphere_TS_lg_erk_lc_n_erk::SWE_Sphere_TS_lg_erk_lc_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		timestepping_order(-1)
{
}



SWE_Sphere_TS_lg_erk_lc_n_erk::~SWE_Sphere_TS_lg_erk_lc_n_erk()
{
}

