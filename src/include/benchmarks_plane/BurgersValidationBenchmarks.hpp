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
#include <sweet/plane/PlaneData.hpp>


class BurgersValidationBenchmarks
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

public:
	static
	double return_u(
			SimulationVariables &i_parameters,
			double x,
			double y
	)
	{
		if (i_parameters.setup.benchmark_scenario_id == 70)
		{
			return sin(2*M_PI*x);
		}

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
			double t=i_parameters.timecontrol.current_simulation_time+i_parameters.timecontrol.current_timestep_size;
			double tmpvar = 0;
			int kmax = i_parameters.sim.f0;
			double eps = 0.1;
			for (int k=1; k<kmax; k++)
			{
				double argument = 2*M_PIl*k*x - M_PIl*k*t + M_PIl*k;
				tmpvar += sin(argument)*eps/sinh(M_PIl*k*eps*0.5);
			}
			tmpvar *= 0.5;
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

	static
	void set_source(
			double i_simulation_time,
			const SimulationVariables &i_parameters,
			bool i_use_staggering,
			PlaneData &io_u_t
	)
	{
		double t = i_simulation_time;
		double tp = 2.0*M_PIl;

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 1/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 57)
		{
			double k = i_parameters.sim.f0;

			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					else
					{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = tp * std::sin(tp*k*x) * std::cos(tp*k*t)
								  + tp * std::sin(tp*k*x) * std::sin(tp*k*t) * std::cos(tp*k*x)*std::sin(tp*k*t)
								  + i_parameters.sim.viscosity * (tp*tp*k * std::sin(tp*k*x) * std::sin(tp*k*t));

					io_data = tmpvar;
				}
			);
		}

		/*
		 * f(t,x,y) = 2*2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)+4*2*PI*sin(2*PI*k*x)*cos(2*PI*k*x)*sin^2(2*PI*k*t)
		 *          - 2*nu(-4*PI^2*k*sin(2*PI*k*x)*sin(2*PI*k*t))
		 * matching to:
		 * u(t,x,y) = 2/k * sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 59)
		{
			double k = i_parameters.sim.f0;


			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = 2*tp * std::sin(tp*k*x) * std::cos(tp*k*t)
								  + 4*tp * std::sin(tp*k*x) * std::sin(tp*k*t) * std::cos(tp*k*x) * std::sin(tp*k*t)
								  + 2*i_parameters.sim.viscosity * (tp*tp*k * std::sin(tp*k*x) * std::sin(tp*k*t));

					io_data = tmpvar;
				}
			);
		}

		/*
		 * f(t,x,y) = 2*PI*sin(2*PI*x)*cos(2*PI*t)+2*PI*sin(2*PI*k*x)*cos(2*PI*k*t)
		 *			+ [sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 *			* [2*PI*cos(2*PI*x)*sin(2*PI*t)+2*PI*cos(2*PI*k*x)*sin(2*PI*k*t)]
		 *          - NU*[-4*PI*PI*sin(2*PI*x)*sin(2*PI*t)
		 *          - 4*PI*PI*k*sin(2*PI*k*x)*sin(2*PI*k*t)]
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)*sin(2*PI*t)+1/k*sin(2*PI*k*x)*sin(2*PI*k*t)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 58)
		{
			double k = i_parameters.sim.f0;

			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = tp * std::sin(tp*x) * std::cos(tp*t) + tp * std::sin(tp*k*x) * std::cos(tp*k*t)
								  + (std::sin(tp*x) * std::sin(tp*t) + 1/k * std::sin(tp*k*x) * std::sin(tp*k*t))
								  * (tp * std::cos(tp*x) * std::sin(tp*t) + tp * std::cos(tp*k*x) * std::sin(tp*k*t))
								  - i_parameters.sim.viscosity * (-tp*tp * std::sin(tp*x) * std::sin(tp*t)
								  - tp*tp*k * std::sin(tp*k*x) * std::sin(tp*k*t));

					io_data = tmpvar;
				}
			);
		}

		/*
		 * f(t,x,y) = 1
		 * matching to:
		 * u(t,x,y) = t
		 */
		if (i_parameters.setup.benchmark_scenario_id == 51)
		{
			io_u_t.physical_set_all(1.0);
		}

		/*
		 * f(t,x,y) = 2*t
		 * matching to:
		 * u(t,x,y) = t^2
		 */
		if (i_parameters.setup.benchmark_scenario_id == 52)
		{
			io_u_t.physical_set_all(2.0*t);
		}

		/*
		 * f(t,x,y) = 3*t^2
		 * matching to:
		 * u(t,x,y) = t^3
		 */
		if (i_parameters.setup.benchmark_scenario_id == 53)
		{
			io_u_t.physical_set_all(3.0*t*t);
		}

		/*
		 * f(t,x,y) = 1000*sin(2*PI*x) + 1000^2*t*sin(2*PI*x)*t*cos(2*PI*x)*2*PI - 1000*NU*(-4*PI*PI*t*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = 1000*t*sin(2*PI*x)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 54)
		{

			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = 1000 * std::sin(tp*x) + 1000*1000*t * std::sin(tp*x) * t * std::cos(tp*x) * tp
								  - 1000 * i_parameters.sim.viscosity * (-tp*tp*t * std::sin(tp*x));

					io_data = tmpvar;
				}
			);
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*t)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*t)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 55)
		{
			io_u_t.physical_set_all(tp*std::cos(tp*t));
		}

		/*
		 * f(t,x,y) = 2*PI*cos(2*PI*k*t)
		 * matching to:
		 * u(t,x,y) = 1/k*sin(2*PI*k*t)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 56)
		{
			double k=i_parameters.sim.f0;
			io_u_t.physical_set_all(tp*std::cos(tp*k*t));
		}

		/*
		 * f(t,x,y) = sin(2*PI*x)*cos(2*PI*x)*2*PI - NU*(-4*PI*PI*sin(2*PI*x))
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 60)
		{

			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = std::sin(tp*x) * std::cos(tp*x) * tp
								  - i_parameters.sim.viscosity * (-tp*tp * std::sin(tp*x));

					io_data = tmpvar;
				}
			);
		}

		/* Test for 2u-grad(u)=F
		 * f(t,x,y) = (1+4*PI*PI)sin(2*PI*x)
		 * matching to:
		 * u(t,x,y) = sin(2*PI*x)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 61)
		{

			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = (1+tp*tp) * std::sin(tp*x);

					io_u_t = tmpvar;
				}
			);
		}

		/*
		 * f(t,x,y) = PI/2* SUM [k*# * (NU*sin(**)*4*PI*k - cos(**) + cos(**)*SUM(sin(**)*#)]
		 * mit: ** = 2*PI*k*x-PI*k*t+PI*k
		 * und: #  = EPS/sinh(0.5*PI*k*EPS)
		 * matching to:
		 * u(t,x,y) = 0.5*SUM_(k=1)^(k_max) sin(2*PI*k*x-PI*k*t+PI*k)*EPS/sinh(0.5*PI*k*EPS)
		 */
		if (i_parameters.setup.benchmark_scenario_id == 62)
		{
			int kmax = i_parameters.sim.f0;
			double eps = 0.1;


			io_u_t.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					// u space
					double x = 0.0;
					if (i_use_staggering)
					{
						x = (((double)i)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}else{
						x = (((double)i+0.5)/(double)i_parameters.disc.res_physical[0])*i_parameters.sim.domain_size[0];
					}
					double tmpvar = 0;
					double A1 = 0;
					double A2 = 0;
					double argument = 0;
					double AA = 0;
					for (int k = 1; k < kmax; k++)
					{
						argument = tp*k*x + M_PIl*k*(1-t);
						AA = eps/sinh(M_PIl*k*eps/2);
						tmpvar += (i_parameters.sim.viscosity*sin(argument)*4*M_PIl*k - cos(argument))*k*AA;
						A1 += cos(argument)*k*AA*M_PIl;
						A2 += sin(argument)*0.5*AA;
					}

					tmpvar *= 0.5*M_PIl;
					tmpvar += A1*A2;
					io_data = tmpvar;
				}
			);
		}

		if (i_parameters.setup.benchmark_scenario_id == 63)
			io_u_t.physical_set_all(0.0);

		if (i_parameters.setup.benchmark_scenario_id == 70)
			io_u_t.physical_set_all(0.0);
	}

	static void printScenarioInformation()
	{
		std::cout << "Available benchmark scenarios:" << std::endl;
		std::cout << "		58 : sin(2*PI*x)*sin(2*PI*t) + sin(2*PI*x*k)*sin(2*PI*t*k)/k" << std::endl;
	}

};

#endif /* SRC_INCLUDE_SWEET_BURGERSVALIDATIONBENCHMARKS_HPP_ */
