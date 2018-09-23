/*
 * SWEPlaneBenchmarks.hpp
 *
 *  Created on: 5 Jul 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SWE_PLANE_BENCHMARKS_HPP_
#define SWE_PLANE_BENCHMARKS_HPP_

#include <iostream>
#include <cassert>
#include <sweet/sweetmath.hpp>
#include <cassert>
#include <stdlib.h>

#include <sweet/SimulationVariables.hpp>


class SWEPlaneBenchmarks_DEPRECATED
{
#if 1
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
			SimulationVariables &i_simVars,
			double x,
			double y,
			int i_variable_id	//< 0:h, 1:u, 2:v, 3:force h, 4:force u, 5:force v, 6:coriolis, 7:divergence, 8: vorticity
	)
	{
		double f = i_simVars.sim.f0;
		double sx = i_simVars.sim.domain_size[0];
		double sy = i_simVars.sim.domain_size[1];

		if (i_simVars.setup.benchmark_scenario_id == 0)
		{
			// radial dam break
			double dx = x-i_simVars.setup.setup_coord_x*sx;
			double dy = y-i_simVars.setup.setup_coord_y*sy;

			if (i_variable_id == 0)	// height
			{
				double radius = i_simVars.setup.radius_scale*sqrt(sx*sx+sy*sy);
				if (dx*dx+dy*dy < radius*radius)
					return i_simVars.sim.h0+1.0;
				else
					return i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 1)
		{
			// Gaussian
			double dx = x-i_simVars.setup.setup_coord_x*sx;
			double dy = y-i_simVars.setup.setup_coord_y*sy;

			if (i_variable_id == 0)
			{
				double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return i_simVars.sim.h0+std::exp(-50.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 2)
		{
			// Steady state (linear and nonlinear) with dominant meridional (y) flow

			if (i_variable_id == 0)
			{
				return std::sin(2.0*M_PI*x/sx) + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					FatalError("f-value is equal to zero!");
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				if (f == 0)
					FatalError("f-value is equal to zero!");
				return i_simVars.sim.gravitation/f*2.0*M_PIl*std::cos(2.0*M_PIl*x/sx)/sx;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 3)
		{
			// Steady state (linear and nonlinear) with dominant zonal (x) flow

			if (i_variable_id == 0)
			{
				return std::sin(2.0*M_PI*y/sy) + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					FatalError("f-value is equal to zero! Cannot run this case scenario.");
				return -i_simVars.sim.gravitation*2.0*M_PI*std::cos(2.0*M_PI*y/sy)/(f*sy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				if (f == 0)
					FatalError("f-value is equal to zero!");
				return 0;

				FatalError("Variable v not available");
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 4)
		{
				double dx = x/sx;
				double dy = y/sy;

			if (i_variable_id == 0)
			{
				return i_simVars.sim.h0 + (std::abs(dx-0.5) < 0.3)*(std::abs(dy-0.5) < 0.1);
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 5)
		{
			double dx = x/sx;
			double dy = y/sy;


			if (i_variable_id == 0)
			{
				return std::sin(6.0*M_PIl*dx)*std::cos(4.0*M_PIl*dy) - (1.0/5.0)*std::cos(4.0*M_PIl*dx)*std::sin(2.0*M_PIl*dy) + i_simVars.sim.h0;
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
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}


		if (i_simVars.setup.benchmark_scenario_id == 6)
		{
			// Gaussian
			double dx = x-i_simVars.setup.setup_coord_x*sx;
			double dy = y-i_simVars.setup.setup_coord_y*sy;

			double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
			double e = 50;
			dx /= radius;
			dy /= radius;


			if (i_variable_id == 0)
			{
				return i_simVars.sim.h0+std::exp(-e*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
			{
				double dh = std::exp(-e*(dx*dx + dy*dy));
				return i_simVars.sim.gravitation/f*e*2.0*dy*dh;
			}
			else if (i_variable_id == 2) // velocity v
			{
				// Gaussian
				double dx = x-i_simVars.setup.setup_coord_x*sx;
				double dy = y-i_simVars.setup.setup_coord_y*sy;

				double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				double e = 50;
				dx /= radius;
				dy /= radius;

				double dh = std::exp(-e*(dx*dx + dy*dy));
				return -i_simVars.sim.gravitation/f*e*2.0*dx*dh;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 8)
		{
			// gaussian in x
			double dx = x-i_simVars.setup.setup_coord_x*sx;


			if (i_variable_id == 0)
			{
				double radius = i_simVars.setup.radius_scale*sx;
				dx /= radius;

				return i_simVars.sim.h0+std::exp(-50.0*(dx*dx));
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 0;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}


		if (i_simVars.setup.benchmark_scenario_id == 9)
		{

			if (i_variable_id == 0)
			{
				return i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				return 1;
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 2;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		if (i_simVars.setup.benchmark_scenario_id == 10)
		{
			// beta plane
			// use e.g. parameters -N 64 -C 0.5 -R 4 -f 0.000001  -g 9.81 -H 1000 -X 100000 -Y 100000 -b 0.0000001 -z -S 1 -s 10

			if (i_variable_id == 0)
			{
				return i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				// Gaussian
				double dx = x-i_simVars.setup.setup_coord_x*sx;
				double dy = y-i_simVars.setup.setup_coord_y*sy;

				double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return 10+std::exp(-50.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		//Forced nonlinear case - trigonometric
		if (i_simVars.setup.benchmark_scenario_id == 13)
		{

			if (i_variable_id == 0)
			{
				return std::cos(2.0*M_PI*x/sx)*std::sin(2.0*M_PI*y/sy) + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				if (f == 0)
					FatalError("f-value is equal to zero! Cannot run this case scenario.");

				double factor = -i_simVars.sim.gravitation*2.0*M_PI/(f*sy);

				return factor*std::cos(2.0*M_PI*x/sx)*std::cos(2.0*M_PI*y/sy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				double factor = -i_simVars.sim.gravitation*2.0*M_PI/(f*sx);
				return factor*std::sin(2.0*M_PI*x/sx)*std::sin(2.0*M_PI*y/sy);
			}
			else if (i_variable_id == 3) // forcing h
			{
				return 0;
			}
			else if (i_variable_id == 4) // forcing u
			{
				double factor = i_simVars.sim.gravitation*2.0*M_PI/(f*sy);
				return -factor*factor*(M_PI/sx)*std::sin(4.0*M_PI*x/sx);
			}
			else if (i_variable_id == 5) // forcing v
			{
				double factor = i_simVars.sim.gravitation*2.0*M_PI/(f*sx);
				return factor*factor*(M_PI/sy)*std::sin(4.0*M_PI*y/sy);
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		//Rotated steady state
		if (i_simVars.setup.benchmark_scenario_id == 14)
		{
			//if ( std::abs(f) < 0.01 )
				//std::cout<<"SWEPlaneBenchmarks Warning: f-value is equal to zero or is very small! Problem may be ill posed. Please be careful or set appropriate -f "<<std::endl;

			double freq = 10.0;

			if ( std::abs(f*sy)+abs(f*sx) < 0.001 )
				FatalError("SWEPlaneBenchmark: f-value*domain is equal to zero or very small! Problem may be ill posed. Please set appropriate -f and -X and -Y");

			if (i_variable_id == 0)
			{
				return std::cos(2.0*M_PI*freq*(x/sx+y/sy)) + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				double factor = i_simVars.sim.gravitation*2.0*M_PI*freq/(f*sy);
				return factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
			}
			else if (i_variable_id == 2) // velocity v
			{
				double factor = -i_simVars.sim.gravitation*2.0*M_PI*freq/(f*sx);
				return factor*std::sin(2.0*M_PI*freq*(x/sx+y/sy));
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}




		// Waves scenario for fixed frequencies
		if (i_simVars.setup.benchmark_scenario_id == 16)
		{
			double param_initial_freq_x_mul=1.0;
			double param_initial_freq_y_mul=1.0;
			double dx = x/i_simVars.sim.domain_size[0]*param_initial_freq_x_mul*M_PIl;
			double dy = y/i_simVars.sim.domain_size[1]*param_initial_freq_y_mul*M_PIl;

			if (i_variable_id == 0) // height
			{
			  return std::sin(2.0*dx)*std::cos(2.0*dy) - (1.0/5.0)*std::cos(2.0*dx)*std::sin(4.0*dy) + i_simVars.sim.h0; // original
			  //return std::sin(2.0*dx)*std::cos(2.0*dy) + i_simVars.sim.h0; //simpler
			  //return std::cos(2.0*M_PI*(x/sx+y/sy)) + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
			  return std::cos(4.0*dx)*std::cos(2.0*dy);//+0.000001*std::cos(100.0*dx)*std::cos(80.0*dy); //simplified
			  //return std::cos(2.0*dx)*std::cos(2.0*dy);//+0.000001*std::cos(40.0*dx)*std::cos(33.0*dy); //simplified
			}
			else if (i_variable_id == 2) // velocity v
			{
			  return std::cos(2.0*dx)*std::cos(4.0*dy); //original
			  //return std::cos(2.0*dx)*std::cos(2.0*dy); //simplified
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else if (i_variable_id == 7) // divergence
			{
				return -2.0*M_PIl*std::sin(2.0*dx+2.0*dy);
			}
			else
			{
				return 0;
			}
		}

		// Diagonal advection
		if (i_simVars.setup.benchmark_scenario_id == 17)
		{
			// Gaussian
			double dx = x-i_simVars.setup.setup_coord_x*sx;
			double dy = y-i_simVars.setup.setup_coord_y*sy;

			if (i_variable_id == 0)
			{
				double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return i_simVars.sim.h0+std::exp(-50.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
			{
				return std::cos(2*M_PI*y);
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 1; //std::sin(2*M_PI*x);
			}
			else if (i_variable_id == 6) // f-term
			{
				return 0; //i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}
		if (i_simVars.setup.benchmark_scenario_id == 18)
		{
			// Preparation for flow over mountain
			double u0 = 20.0; //mean constant velocity
			double g = i_simVars.sim.gravitation; // gravity
			double f = i_simVars.sim.f0; // gravity

			if (i_variable_id == 0)
			{
				return -std::sin(2.0*M_PI*y/sy)*f*u0*(0.5*sy/M_PI)/g + i_simVars.sim.h0;
			}
			else if (i_variable_id == 1) // velocity u
			{
				return u0*std::cos(2.0*M_PI*y/sy);
			}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return f;
			}
			else
			{
				return 0;
			}
		}
		if (i_simVars.setup.benchmark_scenario_id == 19)
		{
			// Rotating Gaussian
			double dx = x-i_simVars.setup.setup_coord_x*sx;
			double dy = y-i_simVars.setup.setup_coord_y*sy;

			if (i_variable_id == 0)
			{
				double radius = i_simVars.setup.radius_scale*sqrt((double)sx*(double)sx+(double)sy*(double)sy);
				dx /= radius;
				dy /= radius;

				return i_simVars.sim.h0+std::exp(-100.0*(dx*dx + dy*dy));
			}
			else if (i_variable_id == 1) // velocity u
					{
				return 10;
					}
			else if (i_variable_id == 2) // velocity v
			{
				return 0;
			}
			else if (i_variable_id == 6) // f-term
			{
				return i_simVars.sim.f0;
			}
			else
			{
				return 0;
			}
		}

		std::cerr << "Invalid setup scenario id " << i_simVars.setup.benchmark_scenario_id << std::endl;
		exit(1);
		return 0;
	}

public:
	static
	double return_h(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 0);
	}

public:
	static
	double return_h_perturbed(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 0) - i_simVars.sim.h0;
	}

	static
	double return_u(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 1);
	}

	static
	double return_v(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 2);
	}

	static
	double return_force_h(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 3);
	}

	static
	double return_force_h_perturbed(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 3) - i_simVars.sim.h0;
	}

	static
	double return_force_u(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 4);
	}

	static
	double return_force_v(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 5);
	}

	static
	double return_f(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 6);
	}

	static
	double return_div(
			SimulationVariables &i_simVars,
			double x,
			double y
	)
	{
		return return_variable_value(i_simVars, x, y, 7);
	}

	//Returns true if initial condition is solution of nonlinear equations, false otherwise
	static
	bool scenario_analytical_init(
				SimulationVariables &i_simVars
	)
	{
		if (i_simVars.setup.benchmark_scenario_id == 0)// radial dam break
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 1) // Gaussian
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 2) // Steady state (linear and nonlinear) with dominant zonal (x) flow
			return true;

		if (i_simVars.setup.benchmark_scenario_id == 3) // Steady state (linear and nonlinear) with dominant meridional (y) flow
			return true;

		if (i_simVars.setup.benchmark_scenario_id == 4)// Square break
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 5) // Trigonometric
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 6) // Gaussian
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 8) // gaussian in x
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 9) //Constant
			return false;

		if (i_simVars.setup.benchmark_scenario_id == 13) //Forced nonlinear case - trigonometric
			return true;

		if (i_simVars.setup.benchmark_scenario_id == 14) //Rotated steady state
			return true;

		if (i_simVars.setup.benchmark_scenario_id == 16) // Waves
			return false;

		return false;
	}

	static void printScenarioInformation()
	{
		std::cout << "Available benchmark scenarios:" << std::endl;
		std::cout << "		0 : radial dam break" << std::endl;
		std::cout << "		1 : Gaussian dam break" << std::endl;
		std::cout << "		2 : balanced state x" << std::endl;
		std::cout << "		3 : balanced state y" << std::endl;
		std::cout << "		4 : square break" << std::endl;
		std::cout << "		5 : trigonometric" << std::endl;
		std::cout << "		6 : Gaussian" << std::endl;
		std::cout << "		8 : Gaussian in x" << std::endl;
		std::cout << "		9 : h=H0, v=0, u=0" << std::endl;
		std::cout << "		13 : Forced nonlinear case - trigonometric" << std::endl;
		std::cout << "		14 : Rotated steady state" << std::endl;
		std::cout << "		16 : Waves" << std::endl;
		std::cout << "		19 : Rotating Gaussian" << std::endl;
	}

};

#endif /* SRC_INCLUDE_SWEET_SWEVALIDATIONBENCHMARKS_HPP_ */
