/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_BENCHMARKS_PLANE_PLANEBENCHMARKSCOMBINED_HPP_
#define SRC_INCLUDE_BENCHMARKS_PLANE_PLANEBENCHMARKSCOMBINED_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/SimulationVariables.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarks.hpp>

class PlaneBenchmarksCombined
{
public:
	static
	void setupInitialConditions(
			PlaneData &o_h,
			PlaneData &o_u,
			PlaneData &o_v,
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
	)
	{
		o_h.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*i_simVars.disc.cell_size[0];
					double y = (double)j*i_simVars.disc.cell_size[1];

					io_data = SWEPlaneBenchmarks::return_h(i_simVars, x, y);
				}
			);

		o_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*i_simVars.disc.cell_size[0];
					double y = (double)j*i_simVars.disc.cell_size[1];

					io_data = SWEPlaneBenchmarks::return_u(i_simVars, x, y);
				}
			);

		o_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (double)i*i_simVars.disc.cell_size[0];
					double y = (double)j*i_simVars.disc.cell_size[1];

					io_data = SWEPlaneBenchmarks::return_v(i_simVars, x, y);
				}
			);
	}

};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
