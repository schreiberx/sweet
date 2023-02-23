/*
 *  Created on: 16 Dec 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef PDE_SWE_PLANE_BENCH_GAUSSIAN_BUMP_HPP__
#define PDE_SWE_PLANE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/plane/Plane.hpp>

#include "../pdeSWEPlane/ShackPDESWEPlaneCoefficients.hpp"
#include <sweet/shacksShared/ShackPlaneDiscretization.hpp>
#include <sweet/shacksShared/ShackSWEPlaneBenchmark.hpp>

/**
 * Setup Gaussian Bump
 *
 * (Formerly -s 1)
 **/
class PDESWEPlaneBenchGaussianBump
{
	PlaneOperators &op;

	ShackPDESWEPlaneCoefficients *shackSimCoeffs;
	ShackPlaneDiscretization *shackDisc;
	ShackSWEPlaneBenchmark *shackBenchmark;
	ShackPDEAdvectionPlaneCoefficients *shackPdeAdvectionPlaneCoeffs;


public:
	PDESWEPlaneBenchGaussianBump(
		sweet::ShackDictionary &io_shackDict,
		PlaneOperators &io_op
	)	:
		op(io_op)
	{
		shackSimCoeffs = io_shackDict.getAutoRegistration<ShackPDESWEPlaneCoefficients>();
		shackDisc = io_shackDict.getAutoRegistration<ShackPlaneDiscretization>();
		shackBenchmark = io_shackDict.getAutoRegistration<ShackSWEPlaneBenchmark>();
		shackPdeAdvectionPlaneCoeffs = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneCoefficients>();
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

		double sx = shackSimCoeffs->plane_domain_size[0];
		double sy = shackSimCoeffs->plane_domain_size[1];

		h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(shackSimCoeffs->plane_domain_size[0]/(double)shackDisc->space_res_physical[0]);
				double y = (double)j*(shackSimCoeffs->plane_domain_size[1]/(double)shackDisc->space_res_physical[1]);

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
				io_data = shackPdeAdvectionPlaneCoeffs->advection_velocity[0];
			}
		);

		v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				io_data = shackPdeAdvectionPlaneCoeffs->advection_velocity[1];
			}
		);

		o_h_pert.loadPlaneDataPhysical(h_pert_phys);
		o_u.loadPlaneDataPhysical(u_phys);
		o_v.loadPlaneDataPhysical(v_phys);

	}
};


#endif
