/*
 * SWE_bench_GaussianBump.hpp
 *
 *  Created on: 16 Dec 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SWE_BENCH_GAUSSIAN_BUMP_HPP__
#define SWE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_Physical.hpp>
#include "PDESWEPlaneBench_BaseInterface.hpp"


/**
 * Setup Gaussian Bump
 */
class PDESWEPlaneBench_GaussianBump	:
		public PDESWEPlaneBench_BaseInterface
{

	double scale;

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

				io_data = std::exp(-50.0*(dx*dx + dy*dy)) * this->scale;
			}
		);

		u_phys.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = 0;
			}
		);

		v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
			io_data = 0;
			}
		);

		o_h_pert.loadPlaneDataPhysical(h_pert_phys);
		o_u.loadPlaneDataPhysical(u_phys);
		o_v.loadPlaneDataPhysical(v_phys);

		return true;
	}

	bool setup(
		sweet::PlaneOperators *io_ops,
		sweet::PlaneData_Config *i_planeDataConfig,
		double i_scale
	)
	{
		PDESWEPlaneBench_BaseInterface::setup(io_ops, i_planeDataConfig);

		this->scale = i_scale;

		return true;
	}


};


#endif
