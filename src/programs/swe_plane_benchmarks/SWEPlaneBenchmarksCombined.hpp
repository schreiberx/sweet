/*
 * SWEPlaneBenchmarksCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP
#define SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP

#include <iostream>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	#include "SWE_bench_Polvani.hpp"
	#include "SWE_bench_MergeVortex.hpp"
	#include "SWE_bench_NormalModes.hpp"
#endif

#include "SWE_bench_UnstableJet.hpp"
#include "SWE_bench_UnstableJetFast.hpp"
#include "SWE_bench_UnstableJetAdv.hpp"
#include "SWE_bench_GaussianBump.hpp"


class SWEPlaneBenchmarksCombined
{
public:
	// Simulation variables
	SimulationVariables *simVars;

	// plane or sphere data config
	const void* ext_forces_data_config;

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
					double sx = simVars->sim.plane_domain_size[0];
					double sy = simVars->sim.plane_domain_size[1];

					// Gaussian
					double dx = i_x-i_center_x*sx;
					double dy = i_y-i_center_y*sy;

					if (dx > 0.5*simVars->sim.plane_domain_size[0])
						dx -= simVars->sim.plane_domain_size[0];
					else if (dx < -0.5*simVars->sim.plane_domain_size[0])
						dx += simVars->sim.plane_domain_size[0];

					if (dy > 0.5*simVars->sim.plane_domain_size[1])
						dy -= simVars->sim.plane_domain_size[1];
					else if (dy < -0.5*simVars->sim.plane_domain_size[1])
						dy += simVars->sim.plane_domain_size[1];

					dx /= sx*simVars->benchmark.object_scale;
					dy /= sy*simVars->benchmark.object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (io_simVars.benchmark.benchmark_name == "")
		{
			SWEETError("SWEPlaneBenchmarksCombined: Benchmark name not given");
		}


#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (io_simVars.benchmark.benchmark_name == "polvani")
		{
			SWE_bench_Polvani swe_polvani(io_simVars, io_op);

			swe_polvani.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (io_simVars.benchmark.benchmark_name == "mergevortex")
		{
			SWE_bench_MergeVortex swe_mergevortex(io_simVars, io_op);

			swe_mergevortex.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (io_simVars.benchmark.benchmark_name == "unstablejet")
		{
			double r = 6.37122e6;
			io_simVars.sim.plane_domain_size[0] = 2.0*M_PI*r;
			io_simVars.sim.plane_domain_size[1] = 2.0*M_PI*r;
			io_simVars.sim.plane_rotating_f0 = 0.00014584;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.h0 = 10000;
			io_op.setup(io_simVars.sim.plane_domain_size, io_simVars.disc.space_use_spectral_basis_diffs);
			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;


			SWE_bench_UnstableJet swe_unstablejet(io_simVars, io_op);

			swe_unstablejet.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (io_simVars.benchmark.benchmark_name == "unstablejet_nobump")
		{
			double r = 6.37122e6;
			io_simVars.sim.plane_domain_size[0] = 2.0*M_PI*r;
			io_simVars.sim.plane_domain_size[1] = 2.0*M_PI*r;
			io_simVars.sim.plane_rotating_f0 = 0.00014584;
			io_simVars.sim.gravitation = 9.80616;
			io_simVars.sim.h0 = 10000;
			io_op.setup(io_simVars.sim.plane_domain_size, io_simVars.disc.space_use_spectral_basis_diffs);
			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;


			SWE_bench_UnstableJet swe_unstablejet(io_simVars, io_op);

			swe_unstablejet.setup(
					o_h_pert,
					o_u,
					o_v,
					false
			);

			return true;
		}
		else if (io_simVars.benchmark.benchmark_name == "unstablejetfast")
			{
				SWE_bench_UnstableJetFast swe_unstablejetfast(io_simVars, io_op);

				swe_unstablejetfast.setup(
						o_h_pert,
						o_u,
						o_v
				);

				return true;
			}

		else if (io_simVars.benchmark.benchmark_name == "unstablejetadv")
		{
			SWE_bench_UnstableJetAdv swe_unstablejetadv(io_simVars, io_op);

			swe_unstablejetadv.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (io_simVars.benchmark.benchmark_name == "normalmodes")
		{
			//PlaneDataConfig *planeDataConfig = o_h_pert.planeDataConfig;

			SWE_bench_NormalModes swe_normalmodes(io_simVars, io_op);

			swe_normalmodes.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#endif
		else if (io_simVars.benchmark.benchmark_name == "gaussian_bump")
		{
			SWE_bench_GaussianBump swe_gaussian_bump(io_simVars, io_op);

			swe_gaussian_bump.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}

		else if (io_simVars.benchmark.benchmark_name == "gaussian_bump_advection")
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
						//*o_plane_data = std::cos(r)*u - std::sin(r)*v;
						*o_plane_data = u*(1.0+std::sin(r));
					}
					else if (i_field_id == 2)
					{
						// v-velocity
						//*o_plane_data = std::sin(r)*u + std::cos(r)*v;
						*o_plane_data = v*(1.0+std::cos(r));
					}

					return;
				}

				SWEETError("Non-existing external field requested!");
				return;
			};

			if (simVars->sim.advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.planeDataConfig;

				// set callback
				io_simVars.benchmark.getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				io_simVars.benchmark.getExternalForcesUserData = this;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, simVars->timecontrol.current_simulation_time, &o_u, simVars->benchmark.getExternalForcesUserData);
				callback_external_forces_advection_field(2, simVars->timecontrol.current_simulation_time, &o_v, simVars->benchmark.getExternalForcesUserData);
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
					double x = (double)i*(simVars->sim.plane_domain_size[0]/(double)simVars->disc.space_res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					io_data = callback_gaussian_bump(center_x, center_y, x, y, exp_fac);
				}
			);
			return true;
		}
		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_0" ||
				io_simVars.benchmark.benchmark_name == "cylinder"
		)
		{
			double sx = simVars->sim.plane_domain_size[0];
			double sy = simVars->sim.plane_domain_size[1];


			o_h_pert.physical_set_zero();
			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.plane_domain_size[0]/(double)io_simVars.disc.space_res_physical[0]);
					double y = (double)j*(io_simVars.sim.plane_domain_size[1]/(double)io_simVars.disc.space_res_physical[1]);

					// radial dam break
					double dx = x-simVars->benchmark.object_coord_x*sx;
					double dy = y-simVars->benchmark.object_coord_y*sy;

					double radius = simVars->benchmark.object_scale*sqrt(sx*sx+sy*sy);
					if (dx*dx+dy*dy < radius*radius)
						io_data = 1.0;
					else
						io_data = 0.0;
				}
			);

			o_u.physical_set_zero();
			o_v.physical_set_zero();

			return true;
		}
		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_1" ||
				io_simVars.benchmark.benchmark_name == "radial_gaussian_bump"
		)
		{
			double sx = simVars->sim.plane_domain_size[0];
			double sy = simVars->sim.plane_domain_size[1];


			o_h_pert.physical_set_zero();
			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(io_simVars.sim.plane_domain_size[0]/(double)io_simVars.disc.space_res_physical[0]);
					double y = (double)j*(io_simVars.sim.plane_domain_size[1]/(double)io_simVars.disc.space_res_physical[1]);

					// radial dam break
					double dx = x-simVars->benchmark.object_coord_x*sx;
					double dy = y-simVars->benchmark.object_coord_y*sy;

					double radius = simVars->benchmark.object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
					dx /= radius;
					dy /= radius;

					io_data = std::exp(-50.0*(dx*dx + dy*dy));
				}
			);

			o_u.physical_set_zero();
			o_v.physical_set_zero();

			return true;
		}
		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_2" ||
				io_simVars.benchmark.benchmark_name == "steady_state_meridional_flow"
		)
		{
			double f = simVars->sim.plane_rotating_f0;
			double sx = simVars->sim.plane_domain_size[0];
			//double sy = simVars->sim.domain_size[1];

			if (io_simVars.sim.plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			o_h_pert.physical_set_zero();
			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)simVars->disc.space_res_physical[0];
					//double y = (double)j*(simVars->sim.domain_size[1]/(double)simVars->disc.res_physical[1]);

					io_data = std::sin(2.0*M_PI*x);
				}
			);

			o_u.physical_set_zero();

			o_v.physical_set_zero();
			o_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)simVars->disc.space_res_physical[0];
					//double y = (double)j*(simVars->sim.domain_size[1]/(double)simVars->disc.res_physical[1]);

					io_data = simVars->sim.gravitation/f*2.0*M_PIl*std::cos(2.0*M_PIl*x)/sx;
				}
			);

			return true;
		}
		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_3" ||
				io_simVars.benchmark.benchmark_name == "steady_state_zonal_flow"
		)
		{
			double f = simVars->sim.plane_rotating_f0;
			//double sx = simVars->sim.domain_size[0];
			double sy = simVars->sim.plane_domain_size[1];

			if (io_simVars.sim.plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			o_h_pert.physical_set_zero();
			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					//double x = (double)i*(simVars->sim.domain_size[0]/(double)simVars->disc.res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					io_data = std::sin(2.0*M_PI*y/sy);
				}
			);

			o_u.physical_set_zero();
			o_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					//double x = (double)i*(simVars->sim.domain_size[0]/(double)simVars->disc.res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					io_data = -simVars->sim.gravitation*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
				}
			);

			o_v.physical_set_zero();

			return true;
		}
		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_4" ||
				io_simVars.benchmark.benchmark_name == "yadda_yadda_whatever_this_is"
		)
		{
			double sx = simVars->sim.plane_domain_size[0];
			double sy = simVars->sim.plane_domain_size[1];

			if (io_simVars.sim.plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			o_h_pert.physical_set_zero();
			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(simVars->sim.plane_domain_size[0]/(double)simVars->disc.space_res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					// radial dam break
					double dx = x-simVars->benchmark.object_coord_x*sx;
					double dy = y-simVars->benchmark.object_coord_y*sy;

					io_data = (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
				}
			);

			o_u.physical_set_zero();
			o_v.physical_set_zero();

			return true;
		}


		else if (
				io_simVars.benchmark.benchmark_name == "benchmark_id_14" ||
				io_simVars.benchmark.benchmark_name == "rotated_steady_state"
		)
		{
			double freq = 10.0;

			double sx = simVars->sim.plane_domain_size[0];
			double sy = simVars->sim.plane_domain_size[1];

			o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(simVars->sim.plane_domain_size[0]/(double)simVars->disc.space_res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					io_data = std::cos(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			o_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(simVars->sim.plane_domain_size[0]/(double)simVars->disc.space_res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					double factor = simVars->sim.gravitation*2.0*M_PI*freq/(simVars->sim.plane_rotating_f0*sy);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			o_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(simVars->sim.plane_domain_size[0]/(double)simVars->disc.space_res_physical[0]);
					double y = (double)j*(simVars->sim.plane_domain_size[1]/(double)simVars->disc.space_res_physical[1]);

					double factor = -simVars->sim.gravitation*2.0*M_PI*freq/(simVars->sim.plane_rotating_f0*sx);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			return true;
		}

		printBenchmarkInformation();
		SWEETError(std::string("Benchmark ")+io_simVars.benchmark.benchmark_name+ " not found (or not available)");


		return false;
	}

	void printBenchmarkInformation()
	{
		std::cout << "Some available benchmark scenarios (--benchmark-name):" << std::endl;
		std::cout << "		polvani : Polvani et al (1994) initial condition" << std::endl;
		std::cout << "		mergevortex : Vortex merging initial conditions from McRae QJ 2014 paper" << std::endl;
		std::cout << "		unstablejet : Unstable Jet test case" << std::endl;
		std::cout << "		gaussian_bump : Gaussian bump" << std::endl;
		std::cout << "		normalmodes : Normal mode initialization" << std::endl;
		std::cout << "Check out more options in src/include/benchmarks_swe_plane/SWEPlaneBenchmarksCombined.hpp" << std::endl;
	}
};



#endif /* SRC_INCLUDE_benchmarks_swe_plane_COMBINED_HPP_ */
