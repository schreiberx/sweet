/*
 *  Created on: 16 Dec 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_ADVECTION_HPP__
#define PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_ADVECTION_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/core/plane/Plane.hpp>

#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneBenchmarks.hpp"
#include "PDEAdvectionPlaneBench_BaseInterface.hpp"

/**
 * Setup Gaussian Bump with advection
 */
class PDEAdvectionPlaneBenchGaussianBumpAdvection	:
		public PDEAdvectionPlaneBench_BaseInterface
{


public:
	bool setupBenchmark(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v
	)
	{

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
};


#endif
