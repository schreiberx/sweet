/*
 * BenchmarkCombined.hpp
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
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
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_h(i_simVars, x, y);
						   }
						   );

	  o_u.physical_update_lambda_array_indices(
						   [&](int i, int j, double &io_data)
						   {
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_u(i_simVars, x, y);
						   }
						   );

	  o_v.physical_update_lambda_array_indices(
						   [&](int i, int j, double &io_data)
						   {
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_v(i_simVars, x, y);
						   }
						   );
	  o_h.physical_update_lambda_array_indices(
						   [&](int i, int j, double &io_data)
						   {
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_h(i_simVars, x, y);
						   }
						   );

	  o_u.physical_update_lambda_array_indices(
						   [&](int i, int j, double &io_data)
						   {
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_u(i_simVars, x, y);
						   }
						   );

	  o_v.physical_update_lambda_array_indices(
						   [&](int i, int j, double &io_data)
						   {
						     double x = (double)i*(i_simVars.sim.domain_size[0]/(double)i_simVars.disc.res_physical[0]);
						     double y = (double)j*(i_simVars.sim.domain_size[1]/(double)i_simVars.disc.res_physical[1]);

						     io_data = SWEPlaneBenchmarks::return_v(i_simVars, x, y);
						   }
						   );

	}

};



#endif /* SRC_INCLUDE_BENCHMARKS_SPHERE_SPHEREBENCHMARKSCOMBINED_HPP_ */
