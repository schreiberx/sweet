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
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "benchmarks/ShackPDESWEPlaneBenchmarks.hpp"
#include "benchmarks/ShackPDESWEPlaneBench_PolvaniBench.hpp"

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	#include "benchmarks/PDESWEPlaneBench_Polvani.hpp"
	#include "benchmarks/PDESWEPlaneBench_MergeVortex.hpp"
	#include "benchmarks/PDESWEPlaneBench_NormalModes.hpp"
#endif

#include "benchmarks/PDESWEPlaneBench_UnstableJet.hpp"
#include "benchmarks/PDESWEPlaneBench_UnstableJetFast.hpp"
#include "benchmarks/PDESWEPlaneBench_UnstableJetAdv.hpp"
#include "benchmarks/PDESWEPlaneBench_GaussianBump.hpp"


#include "ShackPDESWEPlane.hpp"


#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>


class PDESWEPlaneBenchmarksCombined
{
public:
	sweet::ErrorBase error;

	// plane or sphere data config
	const void* ext_forces_data_config;

	sweet::ShackDictionary *shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWEPlane *shackPDESWEPlane;
	ShackPDESWEPlaneBenchmarks *shackPDESWEPlaneBenchmarks;
	ShackPDESWEPlaneBench_PolvaniBench *shackPDESWEPlaneBench_polvaniBenchmark;

	PDESWEPlaneBenchmarksCombined()	:
		shackDict(nullptr),
		shackPlaneDataOps(nullptr),
		shackTimestepControl(nullptr),
		shackPDESWEPlane(nullptr),
		shackPDESWEPlaneBenchmarks(nullptr)
	{
	}


	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	)
	{
		shackDict = io_shackDict;

		shackPlaneDataOps = shackDict->getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = shackDict->getAutoRegistration<sweet::ShackTimestepControl>();
		shackPDESWEPlane = shackDict->getAutoRegistration<ShackPDESWEPlane>();
		shackPDESWEPlaneBenchmarks = shackDict->getAutoRegistration<ShackPDESWEPlaneBenchmarks>();
		shackPDESWEPlaneBench_polvaniBenchmark = shackDict->getAutoRegistration<ShackPDESWEPlaneBench_PolvaniBench>();

		ERROR_CHECK_WITH_RETURN_BOOLEAN(*io_shackDict);

		return true;
	}


	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		return shackRegistration(&io_shackDict);
	}


	bool clear()
	{
		// TODO
		return true;
	}


public:
	bool setupInitialConditions(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v,
			sweet::PlaneOperators *io_ops,				///< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::PlaneDataConfig *io_planeDataConfig
	)
	{
		assert(io_ops != nullptr);
		assert(io_planeDataConfig != nullptr);

		auto callback_gaussian_bump =
				[&](
						double i_center_x, double i_center_y,
						double i_x, double i_y,
						double i_exp_fac
				)
				{
					double sx = shackPlaneDataOps->plane_domain_size[0];
					double sy = shackPlaneDataOps->plane_domain_size[1];

					// Gaussian
					double dx = i_x-i_center_x*sx;
					double dy = i_y-i_center_y*sy;

					if (dx > 0.5*shackPlaneDataOps->plane_domain_size[0])
						dx -= shackPlaneDataOps->plane_domain_size[0];
					else if (dx < -0.5*shackPlaneDataOps->plane_domain_size[0])
						dx += shackPlaneDataOps->plane_domain_size[0];

					if (dy > 0.5*shackPlaneDataOps->plane_domain_size[1])
						dy -= shackPlaneDataOps->plane_domain_size[1];
					else if (dy < -0.5*shackPlaneDataOps->plane_domain_size[1])
						dy += shackPlaneDataOps->plane_domain_size[1];

					dx /= sx*shackPDESWEPlaneBenchmarks->object_scale;
					dy /= sy*shackPDESWEPlaneBenchmarks->object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (shackPDESWEPlaneBenchmarks->benchmark_name == "")
			return error.set("SWEPlaneBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		if (shackPDESWEPlaneBenchmarks->benchmark_name == "polvani")
		{
			PDESWEPlaneBench_Polvani swe_polvani;
			swe_polvani.shackRegistration(shackDict);
			swe_polvani.setup(io_ops, io_planeDataConfig);
			swe_polvani.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "mergevortex")
		{
			PDESWEPlaneBench_MergeVortex swe_mergevortex;
			swe_mergevortex.shackRegistration(shackDict);
			swe_mergevortex.setup(io_ops, io_planeDataConfig);
			swe_mergevortex.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
				);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "unstablejet")
		{
			double r = 6.37122e6;
			shackPlaneDataOps->plane_domain_size[0] = 2.0*M_PI*r;
			shackPlaneDataOps->plane_domain_size[1] = 2.0*M_PI*r;
			shackPDESWEPlane->plane_rotating_f0 = 0.00014584;
			shackPDESWEPlane->gravitation = 9.80616;
			shackPDESWEPlane->h0 = 10000;

			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;
			io_ops->clear();
			io_ops->setup(io_planeDataConfig, shackPlaneDataOps);


			PDESWEPlaneBench_UnstableJet swe_unstablejet;
			swe_unstablejet.shackRegistration(shackDict);
			swe_unstablejet.setup(io_ops, io_planeDataConfig);
			swe_unstablejet.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "unstablejet_nobump")
		{
			double r = 6.37122e6;
			shackPlaneDataOps->plane_domain_size[0] = 2.0*M_PI*r;
			shackPlaneDataOps->plane_domain_size[1] = 2.0*M_PI*r;
			shackPDESWEPlane->plane_rotating_f0 = 0.00014584;
			shackPDESWEPlane->gravitation = 9.80616;
			shackPDESWEPlane->h0 = 10000;
			io_ops->setup(io_planeDataConfig, shackPlaneDataOps);
			std::cout << "WARNING: OVERWRITING SIMULATION PARAMETERS FOR THIS BENCHMARK!" << std::endl;


			PDESWEPlaneBench_UnstableJet swe_unstablejet(false);
			swe_unstablejet.shackRegistration(shackDict);
			swe_unstablejet.setup(io_ops, io_planeDataConfig);
			swe_unstablejet.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "unstablejetfast")
		{
			PDESWEPlaneBench_UnstableJetFast swe_unstablejetfast;
			swe_unstablejetfast.shackRegistration(shackDict);
			swe_unstablejetfast.setup(io_ops, io_planeDataConfig);
			swe_unstablejetfast.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "unstablejetadv")
		{
			PDESWEPlaneBench_UnstableJetAdv swe_unstablejetadv;
			swe_unstablejetadv.shackRegistration(shackDict);
			swe_unstablejetadv.setup(io_ops, io_planeDataConfig);
			swe_unstablejetadv.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "normalmodes")
		{
			PDESWEPlaneBench_NormalModes swe_normalmodes;
			swe_normalmodes.shackRegistration(shackDict);
			swe_normalmodes.setup(io_ops, io_planeDataConfig);
			swe_normalmodes.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#endif
		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "gaussian_bump" || shackPDESWEPlaneBenchmarks->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDESWEPlaneBench_GaussianBump swe_gaussian_bump;
			swe_gaussian_bump.shackRegistration(shackDict);
			swe_gaussian_bump.setup(io_ops, io_planeDataConfig);
			swe_gaussian_bump.setupBenchmark(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}

		else if (shackPDESWEPlaneBenchmarks->benchmark_name == "gaussian_bump_advection")
		{

			auto callback_external_forces_advection_field =
					[](
							int i_field_id,
							double i_simulation_timestamp,
							sweet::PlaneData_Spectral* io_data,			/// planedata or spheredata
							ShackPDESWEPlaneBenchmarks* i_shackBenchmark		/// user data (pointer to this class)
			)
			{
				sweet::PlaneData_Physical plane_data_phys(io_data->planeDataConfig);

				if (i_field_id >= 1 && i_field_id <= 2)
				{
					double u = i_shackBenchmark->advection_velocity[0];
					double v = i_shackBenchmark->advection_velocity[1];

					double r;
					if (i_shackBenchmark->advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/i_shackBenchmark->advection_velocity[2]*2.0*M_PI;

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

					io_data->loadPlaneDataPhysical(plane_data_phys);

					return;
				}

				SWEETError("Non-existing external field requested!");
				return;
			};

			if (shackPDESWEPlaneBenchmarks->advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.planeDataConfig;

				// set callback
				shackPDESWEPlaneBenchmarks->getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				shackPDESWEPlaneBenchmarks->getExternalForcesUserData = shackPDESWEPlaneBenchmarks;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(1, shackTimestepControl->current_simulation_time, &o_u, shackPDESWEPlaneBenchmarks->getExternalForcesUserData);
				callback_external_forces_advection_field(2, shackTimestepControl->current_simulation_time, &o_v, shackPDESWEPlaneBenchmarks->getExternalForcesUserData);
			}
			else
			{
				sweet::PlaneData_Physical u_phys(o_u.planeDataConfig);
				sweet::PlaneData_Physical v_phys(o_u.planeDataConfig);

				u_phys = shackPDESWEPlaneBenchmarks->advection_velocity[0];
				v_phys = shackPDESWEPlaneBenchmarks->advection_velocity[1];

				o_u.loadPlaneDataPhysical(u_phys);
				o_v.loadPlaneDataPhysical(v_phys);
			}

			double center_x = 0.5;
			double center_y = 0.5;
			double exp_fac = 50.0;

			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					io_data = callback_gaussian_bump(center_x, center_y, x, y, exp_fac);
				}
			);
			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			return true;
		}
		else if (
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_0" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "cylinder"
		)
		{
			double sx = shackPlaneDataOps->plane_domain_size[0];
			double sy = shackPlaneDataOps->plane_domain_size[1];


			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWEPlaneBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWEPlaneBenchmarks->object_coord_y*sy;

					double radius = shackPDESWEPlaneBenchmarks->object_scale*sqrt(sx*sx+sy*sy);
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
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_1" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "radial_gaussian_bump"
		)
		{
			double sx = shackPlaneDataOps->plane_domain_size[0];
			double sy = shackPlaneDataOps->plane_domain_size[1];


			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWEPlaneBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWEPlaneBenchmarks->object_coord_y*sy;

					double radius = shackPDESWEPlaneBenchmarks->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
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
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_2" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "steady_state_meridional_flow"
		)
		{
			double f = shackPDESWEPlane->plane_rotating_f0;
			double sx = shackPlaneDataOps->plane_domain_size[0];
			//double sy = shackSim->domain_size[1];

			if (shackPDESWEPlane->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackPlaneDataOps->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = std::sin(2.0*M_PI*x);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			o_u.spectral_set_zero();

			sweet::PlaneData_Physical v_phys(o_v.planeDataConfig);
			v_phys.physical_set_zero();
			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i/(double)shackPlaneDataOps->space_res_physical[0];
					//double y = (double)j*(shackSim->domain_size[1]/(double)disc->res_physical[1]);

					io_data = shackPDESWEPlane->gravitation/f*2.0*M_PI*std::cos(2.0*M_PI*x)/sx;
				}
			);

			o_v.loadPlaneDataPhysical(v_phys);

			return true;
		}
		else if (
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_3" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "steady_state_zonal_flow"
		)
		{
			double f = shackPDESWEPlane->plane_rotating_f0;
			//double sx = shackSim->domain_size[0];
			double sy = shackPlaneDataOps->plane_domain_size[1];

			if (shackPDESWEPlane->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					io_data = std::sin(2.0*M_PI*y/sy);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);

			sweet::PlaneData_Physical u_phys(o_u.planeDataConfig);
			u_phys.physical_set_zero();
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					io_data = -shackPDESWEPlane->gravitation*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
				}
			);

			o_u.loadPlaneDataPhysical(u_phys);

			o_v.spectral_set_zero();

			return true;
		}
		else if (
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_4" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "yadda_yadda_whatever_this_is"
		)
		{
			double sx = shackPlaneDataOps->plane_domain_size[0];
			double sy = shackPlaneDataOps->plane_domain_size[1];

			if (shackPDESWEPlane->plane_rotating_f0 == 0)
				SWEETError("Coriolis = 0!");

			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					// radial dam break
					double dx = x-shackPDESWEPlaneBenchmarks->object_coord_x*sx;
					double dy = y-shackPDESWEPlaneBenchmarks->object_coord_y*sy;

					io_data = (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.spectral_set_zero();
			o_v.spectral_set_zero();

			return true;
		}


		else if (
				shackPDESWEPlaneBenchmarks->benchmark_name == "benchmark_id_14" ||
				shackPDESWEPlaneBenchmarks->benchmark_name == "rotated_steady_state"
		)
		{
			double freq = 10.0;

			double sx = shackPlaneDataOps->plane_domain_size[0];
			double sy = shackPlaneDataOps->plane_domain_size[1];

			sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
			h_pert_phys.physical_set_zero();
			h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					io_data = std::cos(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			sweet::PlaneData_Physical u_phys(o_u.planeDataConfig);
			u_phys.physical_set_zero();
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					double factor = shackPDESWEPlane->gravitation*2.0*M_PI*freq/(shackPDESWEPlane->plane_rotating_f0*sy);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			sweet::PlaneData_Physical v_phys(o_v.planeDataConfig);
			v_phys.physical_set_zero();
			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
					double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

					double factor = -shackPDESWEPlane->gravitation*2.0*M_PI*freq/(shackPDESWEPlane->plane_rotating_f0*sx);
					io_data = factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
				}
			);

			o_h_pert.loadPlaneDataPhysical(h_pert_phys);
			o_u.loadPlaneDataPhysical(u_phys);
			o_v.loadPlaneDataPhysical(v_phys);
			return true;
		}

		printBenchmarkInformation();
		SWEETError(std::string("Benchmark ")+shackPDESWEPlaneBenchmarks->benchmark_name+ " not found (or not available)");


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



#endif
