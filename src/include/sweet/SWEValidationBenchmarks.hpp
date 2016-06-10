/*
 * SWEValidationBenchmarks.hpp
 *
 *  Created on: 5 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_
#define SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_

#include <sweetmath.hpp>

#include <sweet/SimulationVariables.hpp>


class SWEValidationBenchmarks
{
public:
	static
	double return_h(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup.scenario == 0)
		{
			// radial dam break
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			if (dx*dx+dy*dy < radius*radius)
				return i_parameters.setup.h0+1.0;
			else
				return i_parameters.setup.h0;
		}

		if (i_parameters.setup.scenario == 1)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			dx /= radius;
			dy /= radius;

			return i_parameters.setup.h0+std::exp(-50.0*(dx*dx + dy*dy));
		}

		if (i_parameters.setup.scenario == 2)
		{
			return std::sin(2.0*M_PI*x/i_parameters.sim.domain_size[0]) + i_parameters.setup.h0;
		}

		if (i_parameters.setup.scenario == 3)
		{
			return std::sin(2.0*M_PI*y/i_parameters.sim.domain_size[1]) + i_parameters.setup.h0;
		}

		if (i_parameters.setup.scenario == 4)
		{
			double dx = x/i_parameters.sim.domain_size[0];
			double dy = y/i_parameters.sim.domain_size[1];
			return i_parameters.setup.h0 + (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
		}

		if (i_parameters.setup.scenario == 5)
		{
			double dx = x/i_parameters.sim.domain_size[0];
			double dy = y/i_parameters.sim.domain_size[1];

			return std::sin(6.0*M_PIl*dx)*std::cos(4.0*M_PIl*dy) - (1.0/5.0)*std::cos(4.0*M_PIl*dx)*std::sin(2.0*M_PIl*dy) + i_parameters.setup.h0;
		}


		if (i_parameters.setup.scenario == 6)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			double e = 50;
			dx /= radius;
			dy /= radius;

			return i_parameters.setup.h0+std::exp(-e*(dx*dx + dy*dy));
		}

		if (i_parameters.setup.scenario == 8)
		{
			// gaussian in x
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];

			double radius = i_parameters.setup.radius_scale*i_parameters.sim.domain_size[0];
			dx /= radius;

			return i_parameters.setup.h0+std::exp(-50.0*(dx*dx));
		}


		if (i_parameters.setup.scenario == 9)
		{
			return i_parameters.setup.h0;
		}

		if (i_parameters.setup.scenario == 10)
		{
			// beta plane
			// use e.g. parameters -N 64 -C 0.5 -R 4 -f 0.000001  -g 9.81 -H 1000 -X 100000 -Y 100000 -b 0.0000001 -z -S 1 -s 10
			return i_parameters.setup.h0;
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.scenario << std::endl;
		return 0;
	}



	static
	double return_u(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup.scenario == 0)
			return 0;

		if (i_parameters.setup.scenario == 1)
			return 0;

		if (i_parameters.setup.scenario == 2)
		{
			if (i_parameters.sim.f0 == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return 0;
		}

		if (i_parameters.setup.scenario == 3)
		{
			if (i_parameters.sim.f0 == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return -i_parameters.sim.g*2.0*M_PI*std::cos(2.0*M_PI*y/i_parameters.sim.domain_size[1])/(i_parameters.sim.f0*i_parameters.sim.domain_size[1]);
		}

		if (i_parameters.setup.scenario == 4)
		{
			return 0;
		}

		if (i_parameters.setup.scenario == 5)
		{
			double dx = x/i_parameters.sim.domain_size[0];
			double dy = y/i_parameters.sim.domain_size[1];

			return std::cos(6.0*M_PIl*dx)*std::cos(4.0*M_PIl*dy)-4.0*std::sin(6.0*M_PIl*dx)*std::sin(4.0*M_PIl*dy);
		}


		if (i_parameters.setup.scenario == 6)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			double e = 50;
			dx /= radius;
			dy /= radius;

			double dh = std::exp(-e*(dx*dx + dy*dy));
			return i_parameters.sim.g/i_parameters.sim.f0*e*2.0*dy*dh;
		}

		if (i_parameters.setup.scenario == 8)
		{
			return 0;
		}


		if (i_parameters.setup.scenario == 9)
		{
			return 1;
		}


		if (i_parameters.setup.scenario == 10)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			dx /= radius;
			dy /= radius;

			return 10+std::exp(-50.0*(dx*dx + dy*dy));
		}

		if (i_parameters.setup.scenario == 51)
		{
			return i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.scenario == 52)
		{
			return i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.scenario == 53)
		{
			return i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time*i_parameters.timecontrol.current_simulation_time;
		}

		if (i_parameters.setup.scenario == 54)
		{
			return 1000*i_parameters.timecontrol.current_simulation_time*std::sin(2*M_PI*x);
		}

		if (i_parameters.setup.scenario == 55)
		{
			return std::sin(2*M_PI*i_parameters.timecontrol.current_simulation_time);
		}

		if (i_parameters.setup.scenario == 56)
		{
			return std::sin(2*M_PI*i_parameters.timecontrol.current_simulation_time*i_parameters.sim.f0)/i_parameters.sim.f0;
		}

		if (i_parameters.setup.scenario == 57)
		{
			double k=i_parameters.sim.f0;
			double t=i_parameters.timecontrol.current_simulation_time;
			return std::sin(2*M_PI*x*k)*std::sin(2*M_PI*t*k)/k;
		}

		if (i_parameters.setup.scenario == 58)
		{
			double k=i_parameters.sim.f0;
			double t=i_parameters.timecontrol.current_simulation_time;
			return std::sin(2*M_PI*x)*std::sin(2*M_PI*t) + std::sin(2*M_PI*x*k)*std::sin(2*M_PI*t*k)/k;
		}

		if (i_parameters.setup.scenario >= 59 && i_parameters.setup.scenario <= 61)
		{
			return 0;
		}

		if (i_parameters.setup.scenario == 62)
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

		if (i_parameters.setup.scenario == 63)
		{
			double tmpvar = 0;
			tmpvar = sin(2*M_PIl*x)*sin(2*M_PIl*y);
			return tmpvar;
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.scenario << std::endl;
		return 0;
	}



	static
	double return_v(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup.scenario == 0)
			return 0;

		if (i_parameters.setup.scenario == 1)
			return 0;

		if (i_parameters.setup.scenario == 2)
		{
			if (i_parameters.sim.f0 == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return i_parameters.sim.g/i_parameters.sim.f0*2.0*M_PIl*std::cos(2.0*M_PIl*x/i_parameters.sim.domain_size[0])/i_parameters.sim.domain_size[0];
		}

		if (i_parameters.setup.scenario == 3)
		{
			if (i_parameters.sim.f0 == 0)
			{
				std::cerr << "f-value is equal to zero!" << std::endl;
				exit(-1);
			}
			return 0;
		}

		if (i_parameters.setup.scenario == 4)
			return 0;

		if (i_parameters.setup.scenario == 5)
		{
			double dx = x/i_parameters.sim.domain_size[0];
			double dy = y/i_parameters.sim.domain_size[1];

			return std::cos(6.0*M_PIl*dx)*std::cos(6.0*M_PIl*dy);
		}


		if (i_parameters.setup.scenario == 6)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*i_parameters.sim.domain_size[0];
			double dy = y-i_parameters.setup.setup_coord_y*i_parameters.sim.domain_size[1];

			double radius = i_parameters.setup.radius_scale*sqrt((double)i_parameters.sim.domain_size[0]*(double)i_parameters.sim.domain_size[0]+(double)i_parameters.sim.domain_size[1]*(double)i_parameters.sim.domain_size[1]);
			double e = 50;
			dx /= radius;
			dy /= radius;

			double dh = std::exp(-e*(dx*dx + dy*dy));
			return -i_parameters.sim.g/i_parameters.sim.f0*e*2.0*dx*dh;
		}

		if (i_parameters.setup.scenario == 8)
		{
			return 0;
		}

		if (i_parameters.setup.scenario == 9)
		{
			return 2;
		}

		if (i_parameters.setup.scenario == 10)
		{
			return 0;
		}

		if (i_parameters.setup.scenario >= 51 && i_parameters.setup.scenario <= 70)
		{
			return 0;
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.scenario << std::endl;
		return 0;
	}
};

#endif /* SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_ */
