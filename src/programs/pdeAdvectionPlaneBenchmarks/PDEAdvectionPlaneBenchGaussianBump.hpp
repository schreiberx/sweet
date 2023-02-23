/*
 *  Created on: 16 Dec 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__
#define PDE_Advection_PLANE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/plane/Plane.hpp>

#include "../pdeAdvectionPlane/ShackPDEAdvectionPlane.hpp"
#include <sweet/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneBenchmarks.hpp"

/**
 * Setup Gaussian Bump
 *
 * (Formerly -s 1)
 **/
class PDEAdvectionPlaneBenchGaussianBump
{
	PlaneOperators &op;

	ShackPDEAdvectionPlane *shackPDEAdvectionPlane;
	ShackPlaneDataOps *shackDisc;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmark;


public:
	PDEAdvectionPlaneBenchGaussianBump(
		sweet::ShackDictionary &io_shackDict,
		PlaneOperators &io_op
	)	:
		op(io_op)
	{
		shackPDEAdvectionPlane = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlane>();
		shackDisc = io_shackDict.getAutoRegistration<ShackPlaneDataOps>();
		shackBenchmark = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneBenchmarks>();
	}

	void setup(
			PlaneData_Spectral &o_h_pert,
			PlaneData_Spectral &o_u,
			PlaneData_Spectral &o_v
	)
	{

		PlaneData_Physical h_pert_phys(o_h_pert.planeDataConfig);
		PlaneData_Physical u_phys(o_u.planeDataConfig);
		PlaneData_Physical v_phys(o_v.planeDataConfig);

		double sx = shackDisc->plane_domain_size[0];
		double sy = shackDisc->plane_domain_size[1];

		h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(shackDisc->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
				double y = (double)j*(shackDisc->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

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
