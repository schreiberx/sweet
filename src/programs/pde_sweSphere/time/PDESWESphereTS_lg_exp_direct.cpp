/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "PDESWESphereTS_lg_exp_direct.hpp"

#include <iostream>
#include <cassert>
#include <utility>
#include <sweet/expIntegration/REXI.hpp>
#include <sweet/core/sphere/Convert_SphereDataSpectralComplex_to_SphereDataSpectral.hpp>
#include <sweet/core/sphere/Convert_SphereDataSpectral_to_SphereDataSpectralComplex.hpp>
#include <sweet/core/StopwatchBox.hpp>



bool PDESWESphereTS_lg_exp_direct::implementsTimesteppingMethod(
		const std::string &i_timestepping_method
)
{
	timestepping_method = i_timestepping_method;

	timestepping_order = shackPDESWETimeDisc->timestepping_order;
	timestepping_order2 = shackPDESWETimeDisc->timestepping_order2;
	if (i_timestepping_method == "lg_exp_direct")
		return true;

	return false;
}


bool PDESWESphereTS_lg_exp_direct::setup_auto(
		const std::string &i_timestepping_method,
		sweet::SphereOperators *io_ops
)
{
	timestepping_method = i_timestepping_method;

	return setup_main(ops, "phi0");
}


bool PDESWESphereTS_lg_exp_direct::setup_main(
		const sweet::SphereOperators *io_ops,
		const std::string &i_function_name
)
{
	ops = io_ops;
	sphereDataConfig = ops->sphereDataConfig;

	function_name = i_function_name;
	expFunctions.setup(i_function_name);

	return true;
}


std::string PDESWESphereTS_lg_exp_direct::getIDString()
{
	return "l_exp_direct";
}


void PDESWESphereTS_lg_exp_direct::runTimestep(
	const sweet::SphereData_Spectral &i_prog_phi0,
	const sweet::SphereData_Spectral &i_prog_vrt0,
	const sweet::SphereData_Spectral &i_prog_div0,

	sweet::SphereData_Spectral &o_prog_phi0,
	sweet::SphereData_Spectral &o_prog_vrt0,
	sweet::SphereData_Spectral &o_prog_div0,

	double i_fixed_dt,
	double i_simulation_timestamp
)
{
	o_prog_phi0 = i_prog_phi0;
	o_prog_vrt0 = i_prog_vrt0;
	o_prog_div0 = i_prog_div0;

	runTimestep(o_prog_phi0, o_prog_vrt0, o_prog_div0, i_fixed_dt, i_simulation_timestamp);
}



PDESWESphereTS_lg_exp_direct::PDESWESphereTS_lg_exp_direct()
{
}


PDESWESphereTS_lg_exp_direct::~PDESWESphereTS_lg_exp_direct()
{
}


void PDESWESphereTS_lg_exp_direct::runTimestep(
	sweet::SphereData_Spectral &io_prog_phi,
	sweet::SphereData_Spectral &io_prog_vrt,
	sweet::SphereData_Spectral &io_prog_div,

	double i_dt,
	double i_simulation_timestamp
)
{
	assert(shackSphereDataOps != nullptr);

	#if SWEET_BENCHMARK_TIMINGS
		StopwatchBox::getInstance().rexi.start();
		StopwatchBox::getInstance().rexi_timestepping.start();
	#endif


	/*
	 * Using exponential integrators, we must compute an
	 */

	double ir = 1.0/shackSphereDataOps->sphere_radius;
	// avg. geopotential

	double G = -shackPDESWESphere->h0*shackPDESWESphere->gravitation;

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

			l0 = expFunctions.eval(i_dt*(-sqrt_DG))*l0;
			l1 = expFunctions.eval(i_dt*sqrt_DG)*l1;

			phi0 = -G/sqrt_DG * l0 + G/sqrt_DG* l1;
			div0 = l0 + l1;

			idx++;
		}
	}

	#if SWEET_BENCHMARK_TIMINGS
		StopwatchBox::getInstance().rexi_timestepping.stop();
		StopwatchBox::getInstance().rexi.stop();
	#endif
}

