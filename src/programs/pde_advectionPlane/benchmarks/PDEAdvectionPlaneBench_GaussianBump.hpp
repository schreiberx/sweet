/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__
#define PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/core/plane/Plane.hpp>

#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneBenchmarks.hpp"
#include "PDEAdvectionPlaneBench_BaseInterface.hpp"

/**
 * Setup Gaussian Bump
 */
class PDEAdvectionPlaneBenchGaussianBump	:
		public PDEAdvectionPlaneBench_BaseInterface
{

public:
	bool setupBenchmark(
			sweet::PlaneData_Spectral &o_h_pert,
			sweet::PlaneData_Spectral &o_u,
			sweet::PlaneData_Spectral &o_v
	)
	{
		sweet::PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
		sweet::PlaneData_Physical u_phys(o_u.planeDataConfig);
		sweet::PlaneData_Physical v_phys(o_v.planeDataConfig);

		double sx = shackPlaneDataOps->plane_domain_size[0];
		double sy = shackPlaneDataOps->plane_domain_size[1];

		h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(shackPlaneDataOps->plane_domain_size[0]/(double)shackPlaneDataOps->space_res_physical[0]);
				double y = (double)j*(shackPlaneDataOps->plane_domain_size[1]/(double)shackPlaneDataOps->space_res_physical[1]);

				// Gaussian
				double dx = x-shackBenchmarks->object_coord_x*sx;
				double dy = y-shackBenchmarks->object_coord_y*sy;


				double radius = shackBenchmarks->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				io_data = std::exp(-50.0*(dx*dx + dy*dy));
			}
		);

		u_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmarks->advection_velocity[0];
			}
		);

		v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmarks->advection_velocity[1];
			}
		);

		o_h_pert.loadPlaneDataPhysical(h_pert_phys);
		o_u.loadPlaneDataPhysical(u_phys);
		o_v.loadPlaneDataPhysical(v_phys);

		return true;
	}
};


#endif
