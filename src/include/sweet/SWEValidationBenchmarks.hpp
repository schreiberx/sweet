/*
 * SWEValidationBenchmarks.hpp
 *
 *  Created on: 5 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_
#define SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_

#include <sweet/SimulationParameters.hpp>
#include <math.h>


class SWEValidationBenchmarks
{
public:
	static
	double return_h(
			SimulationParameters &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup_scenario == 0)
		{
			// radial dam break
			double dx = x-i_parameters.setup_coord_x;
			double dy = y-i_parameters.setup_coord_y;

			double radius = i_parameters.setup_radius;
			if (dx*dx+dy*dy < radius*radius)
				return i_parameters.setup_h0+1.0;
			else
				return i_parameters.setup_h0;
		}

		if (i_parameters.setup_scenario == 1)
		{
			double dx = x-i_parameters.setup_coord_x;
			double dy = y-i_parameters.setup_coord_y;

			double radius = i_parameters.setup_radius*10;
			dx /= radius;
			dy /= radius;

			return i_parameters.setup_h0+std::exp(-50.0*(dx*dx + dy*dy));
		}

		if (i_parameters.setup_scenario == 2)
		{
			return std::sin(2.0*M_PI*x/i_parameters.sim_domain_length[0]) + i_parameters.setup_h0;
		}

		if (i_parameters.setup_scenario == 3)
		{
			return std::sin(2.0*M_PI*y/i_parameters.sim_domain_length[1])+i_parameters.setup_h0;
		}

		if (i_parameters.setup_scenario == 4)
		{
			double dx = x-i_parameters.setup_coord_x;
			double dy = y-i_parameters.setup_coord_y;

			double radius = i_parameters.setup_radius*10;
			dx /= radius;
			dy /= radius;

			return i_parameters.setup_h0+std::exp(-50.0*(std::max(dx*dx, dy*dy)));
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup_scenario << std::endl;
		return 0;
	}



	static
	double return_u(
			SimulationParameters &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup_scenario == 0)
			return 0;

		if (i_parameters.setup_scenario == 1)
			return 0;

		if (i_parameters.setup_scenario == 2)
		{
			if (i_parameters.sim_f == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return 0;
		}

		if (i_parameters.setup_scenario == 3)
		{
			if (i_parameters.sim_f == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return -2.0*M_PI*std::cos(2.0*M_PI*y/i_parameters.sim_domain_length[1])/(i_parameters.sim_f*i_parameters.sim_domain_length[1]);
		}

		if (i_parameters.setup_scenario == 4)
			return 0;

		std::cerr << "Invalid setup scenario id " << i_parameters.setup_scenario << std::endl;
		return 0;
	}



	static
	double return_v(
			SimulationParameters &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup_scenario == 0)
			return 0;

		if (i_parameters.setup_scenario == 1)
			return 0;

		if (i_parameters.setup_scenario == 2)
		{
			if (i_parameters.sim_f == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return 2.0*M_PI*std::cos(2.0*M_PI*x/i_parameters.sim_domain_length[0])/(i_parameters.sim_f*i_parameters.sim_domain_length[0]);
		}

		if (i_parameters.setup_scenario == 3)
		{
			if (i_parameters.sim_f == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return 0;
		}


		if (i_parameters.setup_scenario == 4)
			return 0;

		std::cerr << "Invalid setup scenario id " << i_parameters.setup_scenario << std::endl;
		return 0;
	}
};

#endif /* SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_ */
