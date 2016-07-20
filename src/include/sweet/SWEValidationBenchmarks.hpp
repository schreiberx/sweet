/*
 * SWEValidationBenchmarks.hpp
 *
 *  Created on: 5 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_
#define SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_

#include <iostream>
#include <cassert>
#include <sweetmath.hpp>

#include <sweet/SimulationVariables.hpp>

#define error(x)	assert(false);

class SWEValidationBenchmarks
{
#if 0
	static
	void error(const char *i_string)
	{
		std::cerr << "ERROR: " << i_string << std::endl;
		assert(false);
		exit(1);
	}
#endif
	static
	double return_variable_value(
			SimulationVariables &i_parameters,
			double x,
			double y,
			int i_variable_id	//< 0:h, 1:u, 2:v, 3:force h, 4:force u, 5:force v
	)
	{
		double f = i_parameters.sim.f0;
		double sx = i_parameters.sim.domain_size[0];
		double sy = i_parameters.sim.domain_size[1];

		if (i_parameters.setup.scenario == 0)
		{
			// radial dam break
			double dx = x-i_parameters.setup.setup_coord_x*sx;
			double dy = y-i_parameters.setup.setup_coord_y*sy;

			if (i_variable_id == 0)	// height
			{
				double radius = i_parameters.setup.radius_scale*sqrt(sx*sx+sy*sy);
				if (dx*dx+dy*dy < radius*radius)
					return i_parameters.setup.h0+1.0;
				else
					return i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 1)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*sx;
			double dy = y-i_parameters.setup.setup_coord_y*sy;

			if (i_variable_id == 0)
			{
				double radius = i_parameters.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return i_parameters.setup.h0+std::exp(-50.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 2)
		{
			// Steady state (linear and nonlinear) with dominant zonal (x) flow

			if (i_variable_id == 0)
			{
				return std::sin(2.0*M_PI*x/sx) + i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					error("f-value is equal to zero!");
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				if (f == 0)
					error("f-value is equal to zero!");
				return i_parameters.sim.g/f*2.0*M_PIl*std::cos(2.0*M_PIl*x/sx)/sx;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 3)
		{
			// Steady state (linear and nonlinear) with dominant meridional (y) flow

			if (i_variable_id == 0)
			{
				return std::sin(2.0*M_PI*y/sy) + i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					error("f-value is equal to zero! Cannot run this case scenario.");
				return -i_parameters.sim.g*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				if (f == 0)
					error("f-value is equal to zero!");
				return 0;

				error("Variable v not available");
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 4)
		{
				double dx = x/sx;
				double dy = y/sy;

			if (i_variable_id == 0)
			{
				return i_parameters.setup.h0 + (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 5)
		{
			double dx = x/sx;
			double dy = y/sy;


			if (i_variable_id == 0)
			{
				return std::sin(6.0*M_PIl*dx)*std::cos(4.0*M_PIl*dy) - (1.0/5.0)*std::cos(4.0*M_PIl*dx)*std::sin(2.0*M_PIl*dy) + i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				double dx = x/sx;
				double dy = y/sy;

				return std::cos(6.0*M_PIl*dx)*std::cos(4.0*M_PIl*dy)-4.0*std::sin(6.0*M_PIl*dx)*std::sin(4.0*M_PIl*dy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				double dx = x/sx;
				double dy = y/sy;

				return std::cos(6.0*M_PIl*dx)*std::cos(6.0*M_PIl*dy);
			}
			else
			{
				return 0;
			}
		}


		if (i_parameters.setup.scenario == 6)
		{
			// Gaussian
			double dx = x-i_parameters.setup.setup_coord_x*sx;
			double dy = y-i_parameters.setup.setup_coord_y*sy;

			double radius = i_parameters.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
			double e = 50;
			dx /= radius;
			dy /= radius;


			if (i_variable_id == 0)
			{
				return i_parameters.setup.h0+std::exp(-e*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
			{
				double dh = std::exp(-e*(dx*dx + dy*dy));
				return i_parameters.sim.g/f*e*2.0*dy*dh;
			}
			else if (i_variable_id == 2) // velocity v
			{
				// Gaussian
				double dx = x-i_parameters.setup.setup_coord_x*sx;
				double dy = y-i_parameters.setup.setup_coord_y*sy;

				double radius = i_parameters.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				double e = 50;
				dx /= radius;
				dy /= radius;

				double dh = std::exp(-e*(dx*dx + dy*dy));
				return -i_parameters.sim.g/f*e*2.0*dx*dh;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 8)
		{
			// gaussian in x
			double dx = x-i_parameters.setup.setup_coord_x*sx;


			if (i_variable_id == 0)
			{
				double radius = i_parameters.setup.radius_scale*sx;
				dx /= radius;

				return i_parameters.setup.h0+std::exp(-50.0*(dx*dx));
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else
			{
				return 0;
			}
		}


		if (i_parameters.setup.scenario == 9)
		{

			if (i_variable_id == 0)
			{
				return i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 1;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 2;
			}
			else
			{
				return 0;
			}
		}

		if (i_parameters.setup.scenario == 10)
		{
			// beta plane
			// use e.g. parameters -N 64 -C 0.5 -R 4 -f 0.000001  -g 9.81 -H 1000 -X 100000 -Y 100000 -b 0.0000001 -z -S 1 -s 10

			if (i_variable_id == 0)
			{
				return i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				// Gaussian
				double dx = x-i_parameters.setup.setup_coord_x*sx;
				double dy = y-i_parameters.setup.setup_coord_y*sy;

				double radius = i_parameters.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return 10+std::exp(-50.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else
			{
				return 0;
			}
		}

		//Forced nonlinear case - trigonometric
		if (i_parameters.setup.scenario == 13)
		{

			if (i_variable_id == 0)
			{
				return std::cos(2.0*M_PI*x/sx)*std::sin(2.0*M_PI*y/sy) + i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					error("f-value is equal to zero! Cannot run this case scenario.");

				double factor = -i_parameters.sim.g*2.0*M_PI/(f*sy);

				return factor*std::cos(2.0*M_PI*x/sx)*std::cos(2.0*M_PI*y/sy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				double factor = -i_parameters.sim.g*2.0*M_PI/(f*sx);
				return factor*std::sin(2.0*M_PI*x/sx)*std::sin(2.0*M_PI*y/sy);
			}
			else if (i_variable_id == 3) // forcing h
			{
				return 0;
			}
			else if (i_variable_id == 4) // forcing u
			{
				double factor = i_parameters.sim.g*2.0*M_PI/(f*sy);
				return -factor*factor*(M_PI/sx)*std::sin(4.0*M_PI*x/sx);
			}
			else if (i_variable_id == 5) // forcing v
			{
				double factor = i_parameters.sim.g*2.0*M_PI/(f*sx);
				return factor*factor*(M_PI/sy)*std::sin(4.0*M_PI*y/sy);
			}
			else
			{
				return 0;
			}
		}

		//Rotated steady state
		if (i_parameters.setup.scenario == 14)
		{
			if (i_variable_id == 0)
			{
				return std::cos(2.0*M_PI*(x/sx+y/sy)) + i_parameters.setup.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					error("f-value is equal to zero! Cannot run this case scenario.");

				double factor = i_parameters.sim.g*2.0*M_PI/(f*sy);
				return factor*std::sin(2.0*M_PI*(x/sx+y/sy));
			}
			else if (i_variable_id == 2) // velocity v
			{
				double factor = -i_parameters.sim.g*2.0*M_PI/(f*sx);
				return factor*std::sin(2.0*M_PI*(x/sx+y/sy));
			}
			else
			{
				return 0;
			}
		}

		std::cerr << "Invalid setup scenario id " << i_parameters.setup.scenario << std::endl;
		exit(1);
		return 0;
	}

public:
	static
	double return_h(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 0);
	}



	static
	double return_u(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 1);
	}



	static
	double return_v(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 2);
	}

	// Forcings for the shallow water equations //
	//------------------------------------------//

	static
	double return_force_h(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 3);
	}

	static
	double return_force_u(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 4);
	}

	static
	double return_force_v(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 5);
	}

	static
	double return_f(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		return return_variable_value(i_parameters, x, y, 6);
	}


	//Returns true is initial condition is solution of nonlinear equations, false otherwise
	static
	bool scenario_analytical_init(
				SimulationVariables &i_parameters
	)
	{
		if (i_parameters.setup.scenario == 0)// radial dam break
			return false;


		if (i_parameters.setup.scenario == 1) // Gaussian
			return false;


		if (i_parameters.setup.scenario == 2) // Steady state (linear and nonlinear) with dominant zonal (x) flow
			return true;

		if (i_parameters.setup.scenario == 3) // Steady state (linear and nonlinear) with dominant meridional (y) flow
			return true;

		if (i_parameters.setup.scenario == 4)// Square break
			return false;


		if (i_parameters.setup.scenario == 5) // Trigonometric
			return false;

		if (i_parameters.setup.scenario == 6) // Gaussian
			return false;

		if (i_parameters.setup.scenario == 8)// gaussian in x
			return false;

		if (i_parameters.setup.scenario == 9)//Constant
			return false;

		if (i_parameters.setup.scenario == 10) // beta plane
			return false;

		if (i_parameters.setup.scenario == 13)//Forced nonlinear case - trigonometric
			return true;

		//Rotated steady state
		if (i_parameters.setup.scenario == 14)
			return true;

		return false;
	}


};

#endif /* SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_ */
