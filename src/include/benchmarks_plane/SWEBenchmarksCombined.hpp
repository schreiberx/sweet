/*
 * PlaneBenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_PLANE_SWEBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_PLANE_SWEBENCHMARKSCOMBINED_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	#include <benchmarks_plane/SWEPolvani.hpp>
#endif



class SWEBenchmarksCombined
{
public:
	static
	bool setupInitialConditions(
			PlaneData &o_h_pert,
			PlaneData &o_u,
			PlaneData &o_v,
			SimulationVariables &io_simVars,	///< Make this IO, since benchmarks can change simulation parameters
			PlaneOperators &io_op				///< Make this IO, since changes in the simulation parameters might require to also update the operators
	)
	{
		if (io_simVars.setup.benchmark_scenario_name == "")
		{
			o_h_pert.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
				double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
				double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

				io_data = SWEPlaneBenchmarks::return_h_perturbed(io_simVars, x, y);
				}
			);

			o_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
					double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

					io_data = SWEPlaneBenchmarks::return_u(io_simVars, x, y);
				}
			);

			o_v.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
					double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

					io_data = SWEPlaneBenchmarks::return_v(io_simVars, x, y);
				}
			);

			return true;
		}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (io_simVars.setup.benchmark_scenario_name == "polvani")
		{
			SWEPolvani swe_polvani(io_simVars, io_op);

			swe_polvani.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
#endif

		printBenchmarkInformation();
		FatalError(std::string("Benchmark ")+io_simVars.setup.benchmark_scenario_name+ " not found (or not availble)");


		return false;
	}

	static void printBenchmarkInformation()
		{
			std::cout << "Available benchmark scenarios (--benchmark):" << std::endl;
			std::cout << "		polvani : Polvani et al (1994) initial condition" << std::endl;
		}
};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
