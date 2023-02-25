/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "SWE_Sphere_TS_l_erk_n_erk.hpp"

#include <sweet/core/SimulationBenchmarkTiming.hpp>


bool SWE_Sphere_TS_l_erk_n_erk::implements_timestepping_method(const std::string &i_timestepping_method
									)
{
	timestepping_method = i_timestepping_method;
	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
	return i_timestepping_method == "l_erk_n_erk";
}

std::string SWE_Sphere_TS_l_erk_n_erk::string_id()
{
	return "l_erk_n_erk";
}


void SWE_Sphere_TS_l_erk_n_erk::run_timestep(
		SphereData_Spectral &io_phi_pert,	///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,	///< prognostic variables

		double i_fixed_dt,
		double i_simulation_timestamp
)
{
	if (timestepping_order == 1)
	{
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order,
				i_simulation_timestamp
			);

		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order2,
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 2)
	{
		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// FULL time step for non-linear part
		timestepping_rk_nonlinear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt,
				timestepping_order2,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);

		// HALF time step for linear part
		timestepping_rk_linear.run_timestep(
				this,
				&SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear,	///< pointer to function to compute euler time step updates
				io_phi_pert, io_vrt, io_div,
				i_fixed_dt*0.5,
				timestepping_order,		/// This must be 2nd order accurate to get overall 2nd order accurate method
				i_simulation_timestamp
			);
	}
	else if (timestepping_order == 0)
	{
		SWEETError("Please specify the timestepping order via --timestepping-order=[int]");
	}
	else
	{
		SWEETError("programs/swe_sphere/SWE_Sphere_TS_l_erk_n_erk.cpp: This order is not yet supported!");
	}
}



void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_linear(
	const SphereData_Spectral &i_phi,	///< prognostic variables
	const SphereData_Spectral &i_vrt,	///< prognostic variables
	const SphereData_Spectral &i_div,	///< prognostic variables

	SphereData_Spectral &o_phi_t,	///< time updates
	SphereData_Spectral &o_vrt_t,	///< time updates
	SphereData_Spectral &o_div_t,	///< time updates

	double i_simulation_timestamp
)
{
	double gh0 = simVars.sim.gravitation * simVars.sim.h0;

	if (!simVars.sim.sphere_use_fsphere)
	{

		/*
		 * See documentation in [sweet]/doc/swe/swe_sphere_formulation/
		 */
		/*
		 * Step 1a
		 */
		sweet::SphereData_Physical ug(i_phi.sphereDataConfig);
		sweet::SphereData_Physical vg(i_phi.sphereDataConfig);
		op.vrtdiv_to_uv(i_vrt, i_div, ug, vg);

		/*
		 * Step 1b
		 */
		sweet::SphereData_Physical tmpg1 = ug*op.fg;
		sweet::SphereData_Physical tmpg2 = vg*op.fg;

		/*
		 * Step 1c
		 */
		op.uv_to_vrtdiv(tmpg1, tmpg2, o_div_t, o_vrt_t);

		/*
		 * Step 1d
		 */
		o_vrt_t *= -1.0;

		/*
		 * Step 1e
		 */
		o_div_t += -op.laplace(i_phi);

		/*
		 * DIV on velocity field
		 */
		o_phi_t = (-gh0)*i_div;
	}
	else
	{
		o_div_t = -op.laplace(i_phi);

		o_vrt_t = -simVars.sim.sphere_fsphere_f0*i_div;
		o_div_t += simVars.sim.sphere_fsphere_f0*i_vrt;

		o_phi_t = (-gh0)*i_div;
	}
}



void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear(
		const SphereData_Spectral &i_phi,	///< prognostic variables
		const SphereData_Spectral &i_vrt,	///< prognostic variables
		const SphereData_Spectral &i_div,	///< prognostic variables

		SphereData_Spectral &o_phi_dt,	///< time updates
		SphereData_Spectral &o_vrt_dt,	///< time updates
		SphereData_Spectral &o_div_dt,	///< time updates

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

	sweet::SphereData_Physical ug(i_phi.sphereDataConfig);
	sweet::SphereData_Physical vg(i_phi.sphereDataConfig);

	sweet::SphereData_Physical vrtg = i_vrt.toPhys();
	sweet::SphereData_Physical divg = i_div.toPhys();
	op.vrtdiv_to_uv(i_vrt, i_div, ug, vg);

	sweet::SphereData_Physical phig = i_phi.toPhys();

	sweet::SphereData_Physical tmpg1 = ug*(vrtg/*+fg*/);
	sweet::SphereData_Physical tmpg2 = vg*(vrtg/*+fg*/);

	op.uv_to_vrtdiv(tmpg1, tmpg2, o_div_dt, o_vrt_dt);

	o_vrt_dt *= -1.0;

	tmpg1 = ug*phig;
	tmpg2 = vg*phig;

	SphereData_Spectral tmpspec(i_phi.sphereDataConfig);
	op.uv_to_vrtdiv(tmpg1,tmpg2, tmpspec, o_phi_dt);

	o_phi_dt *= -1.0;

	sweet::SphereData_Physical tmpg = 0.5*(ug*ug+vg*vg);

	tmpspec = /*phig+*/tmpg;

	o_div_dt += -op.laplace(tmpspec);


#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.stop();
#endif
}



/**
 * This routine is used by other time step implementations
 */
void SWE_Sphere_TS_l_erk_n_erk::euler_timestep_update_nonlinear(
		SphereData_Spectral &io_phi,		///< prognostic variables
		SphereData_Spectral &io_vrt,	///< prognostic variables
		SphereData_Spectral &io_div,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.start();
#endif

	SphereData_Spectral tmp_phi(io_phi.sphereDataConfig);
	SphereData_Spectral tmp_vrt(io_vrt.sphereDataConfig);
	SphereData_Spectral tmp_div(io_div.sphereDataConfig);

	euler_timestep_update_nonlinear(
			io_phi,
			io_vrt,
			io_div,

			tmp_phi,
			tmp_vrt,
			tmp_div,

			i_simulation_timestamp
		);

	io_phi += i_dt*tmp_phi;
	io_vrt += i_dt*tmp_vrt;
	io_div += i_dt*tmp_div;

#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_nonlinearities.stop();
#endif
}






/*
 * Setup
 */
void SWE_Sphere_TS_l_erk_n_erk::setup(
		int i_order,	///< order of RK time stepping method for non-linear parts
		int i_order2	///< order of RK time stepping method for non-linear parts
)
{
	timestepping_order = i_order;
	timestepping_order2 = i_order2;
}



void SWE_Sphere_TS_l_erk_n_erk::setup_auto()
{
	setup(timestepping_order, timestepping_order2);
}



SWE_Sphere_TS_l_erk_n_erk::SWE_Sphere_TS_l_erk_n_erk(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
}



SWE_Sphere_TS_l_erk_n_erk::~SWE_Sphere_TS_l_erk_n_erk()
{
}

