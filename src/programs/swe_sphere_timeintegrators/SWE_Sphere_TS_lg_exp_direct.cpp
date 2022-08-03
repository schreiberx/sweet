/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "SWE_Sphere_TS_lg_exp_direct.hpp"

#include <iostream>
#include <cassert>
#include <utility>
#include <rexi/REXI.hpp>
#include <sweet/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/SimulationBenchmarkTiming.hpp>



bool SWE_Sphere_TS_lg_exp_direct::implements_timestepping_method(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;

	timestepping_order = simVars.disc.timestepping_order;
	timestepping_order2 = simVars.disc.timestepping_order2;
	if (i_timestepping_method == "lg_exp_direct")
		return true;

	if (i_timestepping_method == "l_exp_direct")
		return true;

	return false;
}


std::string SWE_Sphere_TS_lg_exp_direct::string_id()
{
	return "l_exp_direct";
}


void SWE_Sphere_TS_lg_exp_direct::setup_auto()
{
	timestepping_method = simVars.disc.timestepping_method;

	setup("phi0");
}


void SWE_Sphere_TS_lg_exp_direct::run_timestep(
	const SphereData_Spectral &i_prog_phi0,
	const SphereData_Spectral &i_prog_vrt0,
	const SphereData_Spectral &i_prog_div0,

	SphereData_Spectral &o_prog_phi0,
	SphereData_Spectral &o_prog_vrt0,
	SphereData_Spectral &o_prog_div0,

	double i_fixed_dt,
	double i_simulation_timestamp
)
{
	o_prog_phi0 = i_prog_phi0;
	o_prog_vrt0 = i_prog_vrt0;
	o_prog_div0 = i_prog_div0;

	run_timestep(o_prog_phi0, o_prog_vrt0, o_prog_div0, i_fixed_dt, i_simulation_timestamp);
}



SWE_Sphere_TS_lg_exp_direct::SWE_Sphere_TS_lg_exp_direct(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
	simVars(i_simVars),
	simCoeffs(simVars.sim),
	op(i_op),
	sphereDataConfig(i_op.sphereDataConfig)
{
}


SWE_Sphere_TS_lg_exp_direct::~SWE_Sphere_TS_lg_exp_direct()
{
}


/**
 * Setup the REXI
 */
void SWE_Sphere_TS_lg_exp_direct::setup(
		const std::string &i_function_name
)
{
	timestepping_method = "lg_exp";
	function_name = i_function_name;

	rexiFunctions.setup(i_function_name);
}



void SWE_Sphere_TS_lg_exp_direct::run_timestep(
	SphereData_Spectral &io_prog_phi,
	SphereData_Spectral &io_prog_vrt,
	SphereData_Spectral &io_prog_div,

	double i_dt,
	double i_simulation_timestamp
)
{
	#if SWEET_BENCHMARK_TIMINGS
		SimulationBenchmarkTimings::getInstance().rexi.start();
		SimulationBenchmarkTimings::getInstance().rexi_timestepping.start();
	#endif


	/*
	 * Using exponential integrators, we must compute an
	 */

	double ir = 1.0/simVars.sim.sphere_radius;
	// avg. geopotential

	double G = -simCoeffs.h0*simCoeffs.gravitation;

	/*
	 * See doc/rexi/rexi_for_swe_on_nonrotating_sphere.pdf
	 */

	SWEET_THREADING_SPACE_PARALLEL_FOR
	for (int m = 0; m <= sphereDataConfig->spectral_modes_m_max; m++)
	{
		std::size_t idx = sphereDataConfig->getArrayIndexByModes(m, m);
		for (int n = m; n <= sphereDataConfig->spectral_modes_n_max; n++)
		{
			double D = (double)n*((double)n+1.0)*ir*ir;

			if (D == 0)
			{
				idx++;
				continue;
			}

			std::complex<double> &phi0 = io_prog_phi.spectral_space_data[idx];
			std::complex<double> &div0 = io_prog_div.spectral_space_data[idx];

			// result will be imaginary only!
			std::complex<double> sqrt_DG = std::sqrt(std::complex<double>(D*G));

			// Multiply with Q^{-1}
			std::complex<double> l0 = -sqrt_DG/(2*G) * phi0 + 0.5*div0;
			std::complex<double> l1 = +sqrt_DG/(2*G) * phi0 + 0.5*div0;

			l0 = rexiFunctions.eval(i_dt*(-sqrt_DG))*l0;
			l1 = rexiFunctions.eval(i_dt*sqrt_DG)*l1;

			phi0 = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
			div0 = l0 + l1;

			idx++;
		}
	}

	#if SWEET_BENCHMARK_TIMINGS
		SimulationBenchmarkTimings::getInstance().rexi_timestepping.stop();
		SimulationBenchmarkTimings::getInstance().rexi.stop();
	#endif
}

