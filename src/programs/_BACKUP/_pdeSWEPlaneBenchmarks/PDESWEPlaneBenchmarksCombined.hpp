/*
 *  Created on: 30 Nov 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP
#define SRC_INCLUDE_BENCHMARKS_SWE_PLANE_COMBINED_HPP

#include <iostream>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/SimulationVariables.hpp>

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	#include "PDESWEPlaneBenchPolvani.hpp"
	#include "SWE_bench_MergeVortex.hpp"
	#include "SWE_bench_NormalModes.hpp"
#endif

#include "PDESWEPlaneBenchGaussianBump.hpp"

#include "SWE_bench_UnstableJet.hpp"
#include "SWE_bench_UnstableJetFast.hpp"
#include "SWE_bench_UnstableJetAdv.hpp"

#include <sweet/core/shacks/ShackDictionary.hpp>
#include "../pdeSWEPlane/ShackPDESWEPlaneCoefficients.hpp"
#include "../pdeAdvectionPlane/ShackPDEAdvectionPlaneCoefficients.hpp"
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include "../../include/sweet/shacksShared/ShackDiagnostics.hpp"
#include <sweet/core/shacksShared/ShackPlaneDiscretization.hpp>
#include "../../include/sweet/shacksShared/ShackShackMisc.hpp"
#include <sweet/core/shacksShared/ShackSWEPlaneBenchmark.hpp>


class PDESWEPlaneBenchmarksCombined
{
public:
	sweet::ErrorBase error;

	// Simulation variables
	//SimulationVariables *simVars;

	// plane or sphere data config
	const void* ext_forces_data_config;

	sweet::ShackDictionary *shackDict;

	ShackPDESWEPlaneCoefficients *shackSWECoeffs;
	ShackPDEAdvectionPlaneCoefficients *shackAdvectionCoeffs;
	ShackPlaneDiscretization *shackDisc;
	ShackSWEPlaneBenchmark *shackBenchmark;
	ShackTimestepControl *shackTimestepControl;

	PlaneOperators *op;
	sweet::ProgramArguments *programArguments;

	/*
	 * Special function to register shacks for benchmarks.
	 *
	 * This is in particular important for the --help output function to include all benchmarks.
	 */
	void shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		io_shackDict.registerFirstTime<ShackPDESWEBenchPolvani>();
	}


public:
	bool setupInitialConditions(
			PlaneData_Spectral &o_h_pert,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v,
			sweet::ShackDictionary &io_shackDict,
			PlaneOperators &io_op,				///< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::ProgramArguments &io_programArguments
	)
	{
		shackDict = &io_shackDict;

		shackSWECoeffs = io_shackDict.getAutoRegistration<ShackPDESWEPlaneCoefficients>();
		shackAdvectionCoeffs = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneCoefficients>();
		shackDisc = io_shackDict.getAutoRegistration<ShackPlaneDiscretization>();
		shackBenchmark = io_shackDict.getAutoRegistration<ShackSWEPlaneBenchmark>();
		shackTimestepControl = io_shackDict.getAutoRegistration<ShackTimestepControl>();

		if (io_shackDict.error.exists())
			return error.forwardWithPositiveReturn(io_shackDict.error);

		op = &io_op;

		programArguments = &io_programArguments;

		return _setupInitialConditions(
				o_h_pert,
				o_u,
				o_v
			);
	}

public:
	bool _setupInitialConditions(
			PlaneData_Spectral &o_h_pert,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
	)
	{

		auto callback_gaussian_bump =
				[&](
						double i_center_x, double i_center_y,
						double i_x, double i_y,
						double i_exp_fac
				)
				{
					double sx = shackSWECoeffs->plane_domain_size[0];
					double sy = shackSWECoeffs->plane_domain_size[1];

					// Gaussian
					double dx = i_x-i_center_x*sx;
					double dy = i_y-i_center_y*sy;

					if (dx > 0.5*shackSWECoeffs->plane_domain_size[0])
						dx -= shackSWECoeffs->plane_domain_size[0];
					else if (dx < -0.5*shackSWECoeffs->plane_domain_size[0])
						dx += shackSWECoeffs->plane_domain_size[0];

					if (dy > 0.5*shackSWECoeffs->plane_domain_size[1])
						dy -= shackSWECoeffs->plane_domain_size[1];
					else if (dy < -0.5*shackSWECoeffs->plane_domain_size[1])
						dy += shackSWECoeffs->plane_domain_size[1];

					dx /= sx*shackBenchmark->object_scale;
					dy /= sy*shackBenchmark->object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (shackBenchmark->benchmark_name == "")
			return error.set("SWEPlaneBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");


		if (0)
		{
		}
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		else if (shackBenchmark->benchmark_name == "polvani")
		{
			SWEPlaneBenchPolvani swePolvani(*shackDict, *op, *programArguments);

			if (swePolvani.error.exists())
				return error.forwardWithPositiveReturn(swePolvani.error);

			swePolvani.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
#if 0
		else if (shackBenchmark->benchmark_name == "mergevortex")
		{
			SWE_bench_MergeVortex swe_mergevortex(shackDict, *op);

			swe_mergevortex.setup(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (shackBenchmark->benchmark_name == "unstablejet")
		{
			double r = 6.37122e6;
			shackSWECoeffs->plane_domain_size[0] = 2.0*M_PI*r;
			shackSWECoeffs->plane_domain_size[1] = 2.0*M_PI*r;
			shackSWECoeffs->plane_rotating_f0 = 0.00014584;
			shackSWECoeffs->gravitation = 9.80616;
			shackSWECoeffs->h0 = 10000;
			op->setup(shackSWECoeffs->plane_domain_size, shackDisc->space_use_spectral_basis_diffs);
			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;


			SWE_bench_UnstableJet swe_unstablejet(shackDict, *op);

			swe_unstablejet.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackBenchmark->benchmark_name == "unstablejet_nobump")
		{
			double r = 6.37122e6;
			shackSWECoeffs->plane_domain_size[0] = 2.0*M_PI*r;
			shackSWECoeffs->plane_domain_size[1] = 2.0*M_PI*r;
			shackSWECoeffs->plane_rotating_f0 = 0.00014584;
			shackSWECoeffs->gravitation = 9.80616;
			shackSWECoeffs->h0 = 10000;
			op->setup(shackSWECoeffs->plane_domain_size, shackDisc->space_use_spectral_basis_diffs);
			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;


			SWE_bench_UnstableJet swe_unstablejet(shackDict, *op);

			swe_unstablejet.setup(
					o_h_pert,
					o_u,
					o_v,
					false
			);

			return true;
		}
		else if (shackBenchmark->benchmark_name == "unstablejetfast")
			{
				SWE_bench_UnstableJetFast swe_unstablejetfast(shackDict, *op);

				swe_unstablejetfast.setup(
						o_h_pert,
						o_u,
						o_v
				);

				return true;
			}

		else if (shackBenchmark->benchmark_name == "unstablejetadv")
		{
			SWE_bench_UnstableJetAdv swe_unstablejetadv(shackDict, *op);

			swe_unstablejetadv.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackBenchmark->benchmark_name == "normalmodes")
		{
			//PlaneDataConfig *planeDataConfig = o_h_pert.planeDataConfig;

			SWE_bench_NormalModes swe_normalmodes(shackDict, *op);

			swe_normalmodes.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#endif
#endif

#if 1
		else if (shackBenchmark->benchmark_name == "gaussian_bump" || shackBenchmark->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDESWEPlaneBenchGaussianBump swe_gaussian_bump(*shackDict, *op);

			swe_gaussian_bump.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#endif

#if 0
		else if (shackBenchmark->benchmark_name == "gaussian_bump_advection")
		{

			auto callback_external_forces_advection_field =
					[](
							int i_field_id,
							double i_simulation_timestamp,
							void* o_data_void,			/// planedata or spheredata
							void* o_data_user_void		/// user data (pointer to this class)
			)
			{
				PlaneData_Spectral* o_plane_data = (PlaneData_Spectral*)o_data_void;
				PlaneData_Physical plane_data_phys(o_plane_data->planeDataConfig);
				PDESWEPlaneBenchmarksCombined* s = (PDESWEPlaneBenchmarksCombined*)o_data_user_void;

				if (i_field_id >= 1 && i_field_id <= 2)
				{
					double u = s->shackSWECoeffs->advection_velocity[0];
					double v = s->shackSWECoeffs->advection_velocity[1];

					double r;
					if (s->shackSWECoeffs->advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/s->shackSWECoeffs->advection_velocity[2]*2.0*M_PI;

					if (i_field_id == 1)
					{
						// u-velocity
						//*o_plane_data = std::cos(r)*u - std::sin(r)*v;
						plane_data_phys = u*(1.0+std::sin(r));
					}
					else if (i_field_id == 2)
					{
						// v-velocity
						//*o_plane_data = std::sin(r)*u + std::cos(r)*v;
						plane_data_phys = v*(1.0+std::cos(r));
					}

					o_plane_data->loadPlaneDataPhysical(plane_data_phys);

					return;
				}

				SWEETError("Non-existing external field requested!");
				return;
			};

			if (shackSWECoeffs->advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.planeDataConfig;

				// set callback
				shackBenchmark->getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				shackBenchmark->getExternalForcesUserData = this;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, shackTimestepControl->current_simulation_time, &o_u, shackBenchmark->getExternalForcesUserData);
				callback_external_forces_advection_field(2, shackTimestepControl->current_simulation_time, &o_v, shackBenchmark->getExternalForcesUserData);
			}
			else
			{
				PlaneData_Physical u_phys(o_u.planeDataConfig);
				PlaneData_Physical v_phys(o_u.planeDataConfig);

				u_phys = shackSWECoeffs->advection_velocity[0];
				v_phys = shackSWECoeffs->advection_velocity[1];

				o_u.loadPlaneDataPhysical(u_phys);
				o_v.loadPlaneDataPhysical(v_phys);
			}

			double center_x = 0.5;
			double center_y = 0.5;
			double exp_fac = 50.0;

			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					io_data = callback_gaussian_bump(center_x, center_y, x, y, exp_fac);
				}
			);
			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			return true;
		}
		else if (
				shackBenchmark->benchmark_name == "benchmark_id_0" ||
				shackBenchmark->benchmark_name == "cylinder"
		)
		{
			double sx = shackSWECoeffs->plane_domain_size[0];
			double sy = shackSWECoeffs->plane_domain_size[1];


			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					// radial dam break
					double dx = x-shackBenchmark->object_coord_x*sx;
					double dy = y-shackBenchmark->object_coord_y*sy;

					double radius = shackBenchmark->object_scale*sqrt(sx*sx+sy*sy);
					if (dx*dx+dy*dy < radius*radius)
						io_data = 1.0;
					else
						io_data = 0.0;
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.spectral_set_zero();
			o_v.spectral_set_zero();

			return true;
		}
		else if (
				shackBenchmark->benchmark_name == "benchmark_id_1" ||
				shackBenchmark->benchmark_name == "radial_gaussian_bump"
		)
		{
			double sx = shackSWECoeffs->plane_domain_size[0];
			double sy = shackSWECoeffs->plane_domain_size[1];


			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					// radial dam break
					double dx = x-shackBenchmark->object_coord_x*sx;
					double dy = y-shackBenchmark->object_coord_y*sy;

					double radius = shackBenchmark->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
					dx /= radius;
					dy /= radius;

					io_data = std::exp(-50.0*(dx*dx + dy*dy));
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.spectral_set_zero();
			o_v.spectral_set_zero();

			return true;
		}
		else if (
				shackBenchmark->benchmark_name == "benchmark_id_2" ||
				shackBenchmark->benchmark_name == "steady_state_meridional_flow"
		)
		{
			double f = shackSWECoeffs->plane_rotating_f0;
			double sx = shackSim->plane_domain_size[0];
			//double sy = sweCoeffs->domain_size[1];

			if (shackSWECoeffs->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackDisc->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = std::sin(2.0*M_PI*x);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			o_u.spectral_set_zero();

			PlaneData_Physical v_phys(o_v.planeDataConfig);
			v_phys.physical_set_zero();
			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackDisc->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = shackSWECoeffs->gravitation/f*2.0*M_PI*std::cos(2.0*M_PI*x)/sx;
				}
			);

			o_v.loadPlaneDataPhysical(v_phys);

			return true;
		}
		else if (
				shackBenchmark->benchmark_name == "benchmark_id_3" ||
				shackBenchmark->benchmark_name == "steady_state_zonal_flow"
		)
		{
			double f = shackSim->plane_rotating_f0;
			//double sx = sweCoeffs->domain_size[0];
			double sy = shackSWECoeffs->plane_domain_size[1];

			if (shackSWECoeffs->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					//double x = (double)i*(shackSim->domain_size[0]/(double)disc->res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					io_data = std::sin(2.0*M_PI*y/sy);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			PlaneData_Physical u_phys(o_u.planeDataConfig);
			u_phys.physical_set_zero();
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					//double x = (double)i*(shackSim->domain_size[0]/(double)disc->res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					io_data = -shackSWECoeffs->gravitation*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
				}
			);

			o_u.loadPlaneDataPhysical(u_phys);

			o_v.spectral_set_zero();

			return true;
		}
#endif
		else if (
				shackBenchmark->benchmark_name == "benchmark_id_4" ||
				shackBenchmark->benchmark_name == "yadda_yadda_whatever_this_is"
		)
		{
			double sx = shackSWECoeffs->plane_domain_size[0];
			double sy = shackSWECoeffs->plane_domain_size[1];

			if (shackSWECoeffs->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					// radial dam break
					double dx = x-shackBenchmark->object_coord_x*sx;
					double dy = y-shackBenchmark->object_coord_y*sy;

					io_data = (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.spectral_set_zero();
			o_v.spectral_set_zero();

			return true;
		}


		else if (
				shackBenchmark->benchmark_name == "benchmark_id_14" ||
				shackBenchmark->benchmark_name == "rotated_steady_state"
		)
		{
			double freq = 10.0;

			double sx = shackSWECoeffs->plane_domain_size[0];
			double sy = shackSWECoeffs->plane_domain_size[1];

			PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					io_data = std::cos(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			PlaneData_Physical u_phys(o_u.planeDataConfig);
			u_phys.physical_set_zero();
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					double factor = shackSWECoeffs->gravitation*2.0*M_PI*freq/(shackSWECoeffs->plane_rotating_f0*sy);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			PlaneData_Physical v_phys(o_v.planeDataConfig);
			v_phys.physical_set_zero();
			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackSWECoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
					double y = (double)j*(shackSWECoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

					double factor = -shackSWECoeffs->gravitation*2.0*M_PI*freq/(shackSWECoeffs->plane_rotating_f0*sx);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.loadPlaneDataPhysical(u_phys);
			o_v.loadPlaneDataPhysical(v_phys);
			return true;
		}

		printBenchmarkInformation();

		return error.set(std::string("Benchmark ")+shackBenchmark->benchmark_name+ " not found (or not available)");
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



#endif
