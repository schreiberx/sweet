/*
 * ODE_Scalar_TS_interface.hpp
 *
 *  Created on: 08 Jun 2022
 *      Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_

const std::complex<double> I(0.0,1.0);

template <typename T>
class ODE_Scalar_TS_interface
{
private:
	T u_prev;
	std::string model;

protected:
	double param_function_L;
	double param_function_N;

protected:
	T function_L(
			T &i_u,
			double i_dt,
			double i_sim_timestamp
		)
	{
		if (model == "ode1")
			return std::sin(i_u);
		else if (model == "cox_matthews_decay")
			return i_u;
#if SWEET_SCALAR_COMPLEX
		else if (model == "cox_matthews_oscillation")
			return I * i_u;
#endif
		else
			SWEETError("Unknown model " + this->model);
	}

	T function_N(
			T &i_u,
			double i_dt,
			double i_sim_timestamp
		)
	{
		if (model == "ode1")
			return std::sin(i_sim_timestamp);
		else if (model == "cox_matthews_decay")
			return std::sin(i_sim_timestamp);
#if SWEET_SCALAR_COMPLEX
		else if (model == "cox_matthews_oscillation")
			return std::exp(I * i_sim_timestamp);
#endif
		else
			SWEETError("Unknown model " + this->model);
	}

public:
	virtual void run_timestep(
			T &io_u,			///< prognostic variables

			double i_dt,		///< time step size
			double i_sim_timestamp
	) = 0;
	////{
	////	double a = this->param_function_L;
	////	double b = this->param_function_N;

	////	io_y += i_dt * (a * std::sin(io_y) + b * std::sin(i_sim_timestamp));
	////}

#if (SWEET_PARAREAL && SWEET_PARAREAL_SCALAR) || (SWEET_XBRAID && SWEET_XBRAID_SCALAR)
	void run_timestep(
			Parareal_GenericData* io_data,

			double i_dt,		///< time step size
			double i_sim_timestamp
	)
	{
		T u = io_data->get_pointer_to_data_Scalar()->simfields[0];

		run_timestep(
				u,
				i_dt,
				i_sim_timestamp
			);

		io_data->get_pointer_to_data_Scalar()->simfields[0] = u;

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
			double i_b,
			std::string i_model
		)
	{
		this->param_function_L = i_a;
		this->param_function_N = i_b;
		this->model = i_model;
	}

};


#endif
