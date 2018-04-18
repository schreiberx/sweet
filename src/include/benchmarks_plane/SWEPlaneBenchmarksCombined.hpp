/*
 * PlaneBenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_PLANE_SWEPLANEBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_PLANE_SWEPLANEBENCHMARKSCOMBINED_HPP_

#include <iostream>
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



class SWEPlaneBenchmarksCombined
{
public:
	// Simulation variables
	SimulationVariables *simVars;

	// plane or sphere data
	//void* ext_forces_data[3];

	// plane or sphere data config
	const void* ext_forces_data_config;

	// last time stamp to update external forces
	//double ext_forces_last_timestamp[3];



	SWEPlaneBenchmarksCombined()
	{
		//ext_forces_data[0] = nullptr;
		//ext_forces_data[1] = nullptr;
		//ext_forces_data[2] = nullptr;

		//ext_forces_last_timestamp[0] = -1;
		//ext_forces_last_timestamp[1] = -1;
		//ext_forces_last_timestamp[2] = -1;
	}



	~SWEPlaneBenchmarksCombined()
	{
#if 0
		if (ext_forces_data[0] != nullptr)
			delete (PlaneData*)ext_forces_data[0];

		if (ext_forces_data[1] != nullptr)
			delete (PlaneData*)ext_forces_data[1];

		if (ext_forces_data[2] != nullptr)
			delete (PlaneData*)ext_forces_data[2];
#endif
	}




public:
	bool setupInitialConditions(
			PlaneData &o_h_pert,
			PlaneData &o_u,
			PlaneData &o_v,
			SimulationVariables &io_simVars,	///< Make this IO, since benchmarks can change simulation parameters
			PlaneOperators &io_op				///< Make this IO, since changes in the simulation parameters might require to also update the operators
	)
	{
		simVars = &io_simVars;


		auto callback_gaussian_bump =
				[&](
						double i_center_x, double i_center_y,
						double i_x, double i_y,
						double i_exp_fac
				)
				{
					double sx = simVars->sim.domain_size[0];
					double sy = simVars->sim.domain_size[1];

					// Gaussian
					double dx = i_x-i_center_x*sx;
					double dy = i_y-i_center_y*sy;

					if (dx > 0.5*simVars->sim.domain_size[0])
						dx -= simVars->sim.domain_size[0];
					else if (dx < -0.5*simVars->sim.domain_size[0])
						dx += simVars->sim.domain_size[0];

					if (dy > 0.5*simVars->sim.domain_size[1])
						dy -= simVars->sim.domain_size[1];
					else if (dy < -0.5*simVars->sim.domain_size[1])
						dy += simVars->sim.domain_size[1];

					dx /= sx*simVars->setup.radius_scale;
					dy /= sy*simVars->setup.radius_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (io_simVars.setup.benchmark_scenario_name == "")
		{
			std::cout << "WARNING: Using -s [int] is deprecated" << std::endl;
			std::cout << "WARNING: TODO: change to use --benchmark [string] for benchmarks" << std::endl;

			o_h_pert.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
				double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
				double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

				io_data = SWEPlaneBenchmarks_DEPRECATED::return_h_perturbed(io_simVars, x, y);
				}
			);

			o_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
					double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

					io_data = SWEPlaneBenchmarks_DEPRECATED::return_u(io_simVars, x, y);
				}
			);

			o_v.physical_update_lambda_array_indices(
					[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.domain_size[0]/(double)io_simVars.disc.res_physical[0]);
					double y = (double)j*(io_simVars.sim.domain_size[1]/(double)io_simVars.disc.res_physical[1]);

					io_data = SWEPlaneBenchmarks_DEPRECATED::return_v(io_simVars, x, y);
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
#endif
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

		else if (io_simVars.setup.benchmark_scenario_name == "gaussian_bump_advection")
		{

			auto callback_external_forces_advection_field =
					[](
							int i_field_id,
							double i_simulation_timestamp,
							void* o_data_void,			/// planedata or spheredata
							void* o_data_user_void		/// user data (pointer to this class)
			)
			{
				PlaneData* o_plane_data = (PlaneData*)o_data_void;
				SWEPlaneBenchmarksCombined* s = (SWEPlaneBenchmarksCombined*)o_data_user_void;

				if (i_field_id >= 1 && i_field_id <= 2)
				{
//					if (s->ext_forces_data[i_field_id] == nullptr)
//						s->ext_forces_data[i_field_id] = new PlaneData(o_plane_data_config);

//					// set pointer to buffer
//					*o_data_void = s->ext_forces_data[i_field_id];

					// nothing to do?
					//if (s->ext_forces_last_timestamp[i_field_id] == s->simVars->timecontrol.current_simulation_time)
					//	return;

					double u = s->simVars->sim.advection_velocity[0];
					double v = s->simVars->sim.advection_velocity[1];

					double r;
					if (s->simVars->sim.advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/s->simVars->sim.advection_velocity[2]*2.0*M_PI;

					if (i_field_id == 1)
					{
						// u-velocity
						*o_plane_data = std::cos(r)*u - std::sin(r)*v;
					}
					else if (i_field_id == 2)
					{
						// u-velocity
						*o_plane_data = std::sin(r)*u + std::cos(r)*v;
					}

					//s->ext_forces_last_timestamp[i_field_id] = s->simVars->timecontrol.current_simulation_time;
					return;
				}

				FatalError("Non-existing external field requested!");
				return;
			};

			if (simVars->sim.advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.planeDataConfig;

				// external forces time stamp
				//ext_forces_last_timestamp[0] = -1;
				//ext_forces_last_timestamp[1] = -1;
				//ext_forces_last_timestamp[2] = -1;

				// set callback
				io_simVars.sim.getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				io_simVars.sim.getExternalForcesUserData = this;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, simVars->timecontrol.current_simulation_time, &o_u, simVars->sim.getExternalForcesUserData);
				callback_external_forces_advection_field(2, simVars->timecontrol.current_simulation_time, &o_v, simVars->sim.getExternalForcesUserData);
			}
			else
			{
				o_u = simVars->sim.advection_velocity[0];
				o_v = simVars->sim.advection_velocity[1];
			}

			double center_x = 0.5;
			double center_y = 0.5;
			double exp_fac = 50.0;

			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(simVars->sim.domain_size[0]/(double)simVars->disc.res_physical[0]);
					double y = (double)j*(simVars->sim.domain_size[1]/(double)simVars->disc.res_physical[1]);

					io_data = callback_gaussian_bump(center_x, center_y, x, y, exp_fac);
				}
			);
			return true;
		}

		printBenchmarkInformation();
		FatalError(std::string("Benchmark ")+io_simVars.setup.benchmark_scenario_name+ " not found (or not availble)");


		return false;
	}

	void printBenchmarkInformation()
	{
		std::cout << "Available benchmark scenarios (--benchmark):" << std::endl;
		std::cout << "		polvani : Polvani et al (1994) initial condition" << std::endl;
	}
};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
