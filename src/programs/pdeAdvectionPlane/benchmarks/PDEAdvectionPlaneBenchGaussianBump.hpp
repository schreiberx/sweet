/*
 *  Created on: 16 Dec 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__
#define PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/plane/Plane.hpp>

#include <sweet/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneBenchmarks.hpp"

/**
 * Setup Gaussian Bump
 **/
class PDEAdvectionPlaneBenchGaussianBump
{
	sweet::PlaneOperators &op;

	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmark;


public:
	PDEAdvectionPlaneBenchGaussianBump(
		sweet::ShackDictionary &io_shackDict,
		sweet::PlaneOperators &io_op
	)	:
		op(io_op)
	{
		shackPlaneDataOps = io_shackDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackBenchmark = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();
	}

	void setup(
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
				double dx = x-shackBenchmark->object_coord_x*sx;
				double dy = y-shackBenchmark->object_coord_y*sy;


				double radius = shackBenchmark->object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				io_data = std::exp(-50.0*(dx*dx + dy*dy));
			}
		);

		u_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmark->advection_velocity[0];
			}
		);

		v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				io_data = shackBenchmark->advection_velocity[1];
			}
		);

		o_h_pert.loadPlaneDataPhysical(h_pert_phys);
		o_u.loadPlaneDataPhysical(u_phys);
		o_v.loadPlaneDataPhysical(v_phys);

	}
};


#endif
