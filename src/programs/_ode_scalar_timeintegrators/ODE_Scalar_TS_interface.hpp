/*
 * ODE_Scalar_TS_interface.hpp
 *
 *  Created on: 08 Jun 2022
 *      Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_



class ODE_Scalar_TS_interface
{
private:
	double u_prev;
	double param_parareal_function_a;
	double param_parareal_function_b;


public:
	void runTimestep(
			double &io_y,			///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		double a = this->param_parareal_function_a;
		double b = this->param_parareal_function_b;

		io_y += i_dt * (a * std::sin(io_y) + b * std::sin(i_sim_timestamp));
	}

#if (SWEET_PARAREAL && SWEET_PARAREAL_SCALAR) || (SWEET_XBRAID && SWEET_XBRAID_SCALAR)
	void runTimestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		double y = io_data->get_pointer_to_data_Scalar()->simfields[0];

		runTimestep(y,
				i_dt,
				i_sim_timestamp
			);

		io_data->get_pointer_to_data_Scalar()->simfields[0] = y;

	}

	// for parareal SL (not needed here)
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{
		u_prev = i_data->get_pointer_to_data_Scalar()->simfields[0];
	};
#endif


	void setup(
			double i_a,
			double i_b
		)
	{
		this->param_parareal_function_a = i_a;
		this->param_parareal_function_b = i_b;
	}

};


#endif
