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

#include "ShackPDEAdvectionPlaneBenchmarks.hpp"
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include "PDEAdvectionPlaneBenchGaussianBump.hpp"



class PDEAdvectionPlaneBenchmarksCombined
{
public:
	sweet::ErrorBase error;

	// Simulation variables
	//SimulationVariables *simVars;

	// plane or sphere data config
	const void* ext_forces_data_config;

	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmarks;

	sweet::PlaneOperators *op;

	PDEAdvectionPlaneBenchmarksCombined()	:
		ext_forces_data_config(nullptr),
		shackPlaneDataOps(nullptr),
		shackTimestepControl(nullptr),
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
		shackPlaneDataOps = io_shackDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = io_shackDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackBenchmarks = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}

	void clear()
	{
		shackPlaneDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackBenchmarks = nullptr;
	}


public:
	bool setupInitialConditions(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v,
			sweet::PlaneOperators &io_op,				///< Make this IO, since changes in the simulation parameters might require to also update the operators
			sweet::ShackDictionary &io_shackDict
	)
	{
		op = &io_op;

		return _setupInitialConditions(
				o_h_pert,
				o_u,
				o_v,
				io_shackDict
			);
	}

public:
	bool _setupInitialConditions(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v,
			sweet::ShackDictionary &io_shackDict
	)
	{
		if (!shackBenchmarks->validateNonzeroAdvection())
			return error.forwardWithPositiveReturn(shackBenchmarks->error);

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

					dx /= sx*shackBenchmarks->object_scale;
					dy /= sy*shackBenchmarks->object_scale;

					return std::exp(-i_exp_fac*(dx*dx + dy*dy));
				};

		if (shackBenchmarks->benchmark_name == "")
			return error.set("SWEPlaneBenchmarksCombined: Benchmark name not given, use --benchmark-name=[name]");


		if (shackBenchmarks->benchmark_name == "gaussian_bump" || shackBenchmarks->benchmark_name == "gaussian_bump_phi_pint")
		{
			PDEAdvectionPlaneBenchGaussianBump gaussian_bump(io_shackDict, *op);

			gaussian_bump.setup(
					o_h_pert,
					o_u,
					o_v
			);

			return true;
		}
#if 1
		else if (shackBenchmarks->benchmark_name == "gaussian_bump_advection")
		{
			auto callback_external_forces_advection_field =
					[](
							int i_field_id,
							double i_simulation_timestamp,
							sweet::PlaneData_Spectral* o_plane_data,			/// planedata or spheredata
							ShackPDEAdvectionPlaneBenchmarks* o_data_user_void		/// user data (pointer to this class)
			)
			{
				sweet::PlaneData_Physical plane_data_phys(o_plane_data->planeDataConfig);
				ShackPDEAdvectionPlaneBenchmarks* shackBenchmarks = (ShackPDEAdvectionPlaneBenchmarks*)o_data_user_void;

				if (i_field_id >= 1 && i_field_id <= 2)
				{
					double u = shackBenchmarks->advection_velocity[0];
					double v = shackBenchmarks->advection_velocity[1];

					double r;
					if (shackBenchmarks->advection_velocity[2] == 0)
						r = 0;
					else
						r = i_simulation_timestamp/shackBenchmarks->advection_velocity[2]*2.0*M_PI;

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

			if (shackBenchmarks->advection_velocity[2] != 0)
			{
				// backup data config
				ext_forces_data_config = o_h_pert.planeDataConfig;

				// set callback
				shackBenchmarks->getExternalForcesCallback = callback_external_forces_advection_field;

				// set user data to this class
				shackBenchmarks->getExternalForcesUserData = this;

				// setup velocities with initial time stamp
				callback_external_forces_advection_field(
						1,
						shackTimestepControl->current_simulation_time,
						&o_u,
						shackBenchmarks
					);

				callback_external_forces_advection_field(
						2,
						shackTimestepControl->current_simulation_time,
						&o_v,
						shackBenchmarks
					);
			}
			else
			{
				sweet::PlaneData_Physical u_phys(o_u.planeDataConfig);
				sweet::PlaneData_Physical v_phys(o_u.planeDataConfig);

				u_phys = shackBenchmarks->advection_velocity[0];
				v_phys = shackBenchmarks->advection_velocity[1];

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
				shackBenchmarks->benchmark_name == "benchmark_id_0" ||
				shackBenchmarks->benchmark_name == "cylinder"
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
					double dx = x-shackBenchmarks->object_coord_x*sx;
					double dy = y-shackBenchmarks->object_coord_y*sy;

					double radius = shackBenchmarks->object_scale*sqrt(sx*sx+sy*sy);
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
				shackBenchmarks->benchmark_name == "benchmark_id_1" ||
				shackBenchmarks->benchmark_name == "radial_gaussian_bump"
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
					double dx = x-shackBenchmarks->object_coord_x*sx;
					double dy = y-shackBenchmarks->object_coord_y*sy;

					double radius = shackBenchmarks->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
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
#endif

		printBenchmarkInformation();

		return error.set(std::string("Benchmark ")+shackBenchmarks->benchmark_name+ " not found (or not available)");
	}

	void printBenchmarkInformation()
	{
		std::cout << "Some available benchmark scenarios (--benchmark-name):" << std::endl;
		std::cout << "		gaussian_bump : Gaussian bump" << std::endl;
	}
};



#endif /* SRC_INCLUDE_benchmarks_swe_plane_COMBINED_HPP_ */
