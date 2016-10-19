/*
 * BurgersValidationBenchmarks.hpp
 *
 *  Created on: 29 Jun 2016
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */
#ifndef SRC_INCLUDE_SWEET_BURGERSVALIDATIONBENCHMARKS_HPP_
#define SRC_INCLUDE_SWEET_BURGERSVALIDATIONBENCHMARKS_HPP_

#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>


class BurgersValidationBenchmarks
{
public:
	static
	double return_u(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{

		if (i_parameters.setup.benchmark_scenario_id == 51)
		{
			return i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.benchmark_scenario_id == 52)
		{
			return i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.benchmark_scenario_id == 53)
		{
			return i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.benchmark_scenario_id == 54)
		{
			return 1000*i_parameters.timecontrol.current_simulation_time*std::sin(2*M_PI*x);
		}

		if (i_parameters.setup.benchmark_scenario_id == 55)
		{
			return std::sin(2*M_PI*i_parameters.timecontrol.current_simulation_time);
		}

		if (i_parameters.setup.benchmark_scenario_id == 56)
		{
			return std::sin(2*M_PI*i_parameters.timecontrol.current_simulation_time*i_parameters.sim.f0)/i_parameters.sim.f0;
		}

		if (i_parameters.setup.benchmark_scenario_id == 57)
		{
			double k=i_parameters.sim.f0;
			double t=i_parameters.timecontrol.current_simulation_time;
			return std::sin(2*M_PI*x*k)*std::sin(2*M_PI*t*k)/k;
		}

		if (i_parameters.setup.benchmark_scenario_id == 58)
		{
			double k=i_parameters.sim.f0;
			double t=i_parameters.timecontrol.current_simulation_time;
			return std::sin(2*M_PI*x)*std::sin(2*M_PI*t) + std::sin(2*M_PI*x*k)*std::sin(2*M_PI*t*k)/k;
		}

		if (i_parameters.setup.benchmark_scenario_id >= 59 && i_parameters.setup.benchmark_scenario_id <= 61)
		{
			return 0;
		}

		if (i_parameters.setup.benchmark_scenario_id == 62)
		{
			double t=i_parameters.timecontrol.current_simulation_time;
			double tmpvar = 0;
			int kmax = 100;
			double eps = 0.1;
			for (int k=1; k<kmax; k++)
			{
				double argument = 2*M_PIl*k*x - M_PIl*k*t + M_PIl*k;
				tmpvar += sin(argument)*eps/sinh(M_PIl*k*eps*0.5);
			}
			tmpvar *= 0.25;
			return tmpvar;
		}

		if (i_parameters.setup.benchmark_scenario_id == 63)
		{
			double tmpvar = 0;
			tmpvar = sin(2*M_PIl*x)*sin(2*M_PIl*y);
			return tmpvar;
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.benchmark_scenario_id << std::endl;
		exit(1);
		return 0;
	}



	static
	double return_v(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup.benchmark_scenario_id >= 51 && i_parameters.setup.benchmark_scenario_id <= 70)
		{
			return 0;
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.benchmark_scenario_id << std::endl;
		exit(1);
		return 0;
	}

};

#endif /* SRC_INCLUDE_SWEET_BURGERSVALIDATIONBENCHMARKS_HPP_ */
