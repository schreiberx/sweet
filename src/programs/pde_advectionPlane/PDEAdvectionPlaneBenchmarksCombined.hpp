/*
 *  Created on: 30 Nov 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP
#define SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP

#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "benchmarks/ShackPDEAdvectionPlaneBenchmarks.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include "benchmarks/PDEAdvectionPlaneBench_GaussianBump.hpp"
#include "benchmarks/PDEAdvectionPlaneBench_GaussianBumpAdvection.hpp"
#include "benchmarks/PDEAdvectionPlaneBench_Cylinder.hpp"
#include "benchmarks/PDEAdvectionPlaneBench_RadialGaussianBump.hpp"



class PDEAdvectionPlaneBenchmarksCombined
{
public:
	sweet::ErrorBase error;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmarks;
	sweet::PlaneOperators *op;

	PDEAdvectionPlaneBenchmarksCombined()	:
		shackBenchmarks(nullptr),
		op(nullptr)
	{

	}


	/*
	 * Special function to register shacks for benchmarks.
	 *
	 * This is in particular important for the --help output function to include all benchmarks.
	 */
	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackBenchmarks = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}

	void clear()
	{
		shackBenchmarks = nullptr;
	}


public:

public:
	bool setupInitialConditions(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v,
			sweet::PlaneOperators &io_op,				///< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::ShackDictionary &io_shackDict
	)
	{
		if (!shackBenchmarks->validateNonzeroAdvection())
			return error.forwardWithPositiveReturn(shackBenchmarks->error);

		if (shackBenchmarks->benchmark_name == "")
			return error.set("SWEPlaneBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");


		if (shackBenchmarks->benchmark_name == "gaussian_bump" || shackBenchmarks->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDEAdvectionPlaneBenchGaussianBump gaussian_bump;

			gaussian_bump.shackRegistration(&io_shackDict);

			gaussian_bump.setup(op);

			gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackBenchmarks->benchmark_name == "gaussian_bump_advection")
		{
			PDEAdvectionPlaneBenchGaussianBumpAdvection gaussian_bump_adv;

			gaussian_bump_adv.shackRegistration(&io_shackDict);

			gaussian_bump_adv.setup(op);

			gaussian_bump_adv.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (
				shackBenchmarks->benchmark_name == "benchmark_id_0" ||
				shackBenchmarks->benchmark_name == "cylinder"
		)
		{

			PDEAdvectionPlaneBenchCylinder gaussian_bump_cylinder;

			gaussian_bump_cylinder.shackRegistration(&io_shackDict);

			gaussian_bump_cylinder.setup(op);

			gaussian_bump_cylinder.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (
				shackBenchmarks->benchmark_name == "benchmark_id_1" ||
				shackBenchmarks->benchmark_name == "radial_gaussian_bump"
		)
		{

			PDEAdvectionPlaneBenchRadialGaussianBump gaussian_radial_gaussian_bump;

			gaussian_radial_gaussian_bump.shackRegistration(&io_shackDict);

			gaussian_radial_gaussian_bump.setup(op);

			gaussian_radial_gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}

		shackBenchmarks->printProgramArguments();

		return error.set(std::string("Benchmark '")+shackBenchmarks->benchmark_name+ "' not found (or not available)");
	}
};



#endif
