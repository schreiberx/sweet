/*
 * SWE_bench_GaussianBump.hpp
 *
 *  Created on: 16 Dec 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */
#ifndef SWE_BENCH_GAUSSIAN_BUMP_HPP__
#define SWE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <cmath>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>


/**
 * Setup Gaussian Bump
 *
 * (Formerly -s 1)
 **/
class SWE_bench_GaussianBump
{
	SimulationVariables &simVars;

	PlaneOperators &op;



public:
	SWE_bench_GaussianBump(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
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

		double sx = simVars.sim.plane_domain_size[0];
		double sy = simVars.sim.plane_domain_size[1];

		h_pert_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(simVars.sim.plane_domain_size[0]/(double)simVars.disc.space_res_physical[0]);
				double y = (double)j*(simVars.sim.plane_domain_size[1]/(double)simVars.disc.space_res_physical[1]);

				// Gaussian
				double dx = x-simVars.benchmark.object_coord_x*sx;
				double dy = y-simVars.benchmark.object_coord_y*sy;


				double radius = simVars.benchmark.object_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				io_data = std::exp(-50.0*(dx*dx + dy*dy));
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

	}
};


#endif
