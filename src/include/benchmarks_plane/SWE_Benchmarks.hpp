/*
 * PlaneBenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_PLANE_SWE_BENCHMARKS_HPP_
#define SRC_INCLUDE_BENCHMARKS_PLANE_SWE_BENCHMARKS_HPP_

#include <benchmarks_plane/SWE_bench_PlaneBenchmarks_DEPRECATED.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	#include <benchmarks_plane/SWE_bench_Polvani.hpp>
	#include <benchmarks_plane/SWE_bench_MergeVortex.hpp>
#endif
#include <benchmarks_plane/SWE_bench_UnstableJet.hpp>
#include <benchmarks_plane/SWE_bench_UnstableJetFast.hpp>
#include <benchmarks_plane/SWE_bench_UnstableJetAdv.hpp>
#include <benchmarks_plane/SWE_bench_GaussianBump.hpp>


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
			std::cout << "WARNING: Using -s [int] is deprecated" << std::endl;
			std::cout << "WARNING: TODO: change to use --benchmark [string] for benchmarks" << std::endl;

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
			SWE_bench_Polvani swe_polvani(io_simVars, io_op);

			swe_polvani.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (io_simVars.setup.benchmark_scenario_name == "mergevortex")
		{
			SWE_bench_MergeVortex swe_mergevortex(io_simVars, io_op);

			swe_mergevortex.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (io_simVars.setup.benchmark_scenario_name == "unstablejet")
		{
			SWE_bench_UnstableJet swe_unstablejet(io_simVars, io_op);

			swe_unstablejet.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (io_simVars.setup.benchmark_scenario_name == "unstablejetfast")
			{
				SWE_bench_UnstableJetFast swe_unstablejetfast(io_simVars, io_op);

				swe_unstablejetfast.setup(
						o_h_pert,
						o_u,
						o_v
				);

				return true;
			}

		else if (io_simVars.setup.benchmark_scenario_name == "unstablejetadv")
		{
			SWE_bench_UnstableJetAdv swe_unstablejetadv(io_simVars, io_op);

			swe_unstablejetadv.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}


		else if (io_simVars.setup.benchmark_scenario_name == "gaussian_bump")
		{
			SWE_bench_GaussianBump swe_gaussian_bump(io_simVars, io_op);

			swe_gaussian_bump.setup(
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
