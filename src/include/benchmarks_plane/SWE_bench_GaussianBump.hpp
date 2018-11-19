/*
 * SWE_bench_GaussianBump.hpp
 *
 *  Created on: 16 Dec 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SWE_BENCH_GAUSSIAN_BUMP_HPP__
#define SWE_BENCH_GAUSSIAN_BUMP_HPP__


#include <stdlib.h>
#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>


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
			PlaneData &o_h_pert,
			PlaneData &o_u,
			PlaneData &o_v
	)
	{
		double sx = simVars.sim.domain_size[0];
		double sy = simVars.sim.domain_size[1];

		o_h_pert.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
				double x = (double)i*(simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0]);
				double y = (double)j*(simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1]);

				// Gaussian
				double dx = x-simVars.benchmark.initial_condition_setup_coord_x*sx;
				double dy = y-simVars.benchmark.initial_condition_setup_coord_y*sy;


				double radius = simVars.benchmark.initial_condition_radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				io_data = std::exp(-50.0*(dx*dx + dy*dy));
			}
		);

		o_u.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = 0;
			}
		);

		o_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
			{
			io_data = 0;
			}
		);
	}
};


#endif
