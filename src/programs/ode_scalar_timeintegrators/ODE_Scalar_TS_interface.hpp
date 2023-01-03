/*
 * ODE_Scalar_TS_interface.hpp
 *
 *  Created on: 08 Jun 2022
 *      Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TS_INTERFACE_HPP_

////////////////////////////////#if (!SWEET_PARAREAL) && (!SWEET_XBRAID)
////////////////////////////////const std::complex<double> I(0.0,1.0);
////////////////////////////////#endif

#include <cassert>

template <typename T>
class ODE_Scalar_TS_interface
{
protected:
	ScalarDataArray u_prev;
	std::string model;

protected:
	std::size_t N;
	std::vector<double> param_function_L;
	std::vector<double> param_function_N;
	std::vector<double> param_function_extra;

protected:
	ScalarDataArray lambda_L(
			double i_dt,
			double i_sim_timestamp
		)
	{
		ScalarDataArray out;

		out.setup(this->N);
		if (this->model == "ode1")
			for (std::size_t i = 0; i < this->N; i++)
				out[i] = std::sin(i_sim_timestamp) * this->param_function_L[i];
		else if (this->model == "cox_matthews_decay")
			for (std::size_t i = 0; i < this->N; i++)
				out[i] = this->param_function_L[i];
#if SWEET_SCALAR_COMPLEX
		else if (this->model == "cox_matthews_oscillation")
			for (std::size_t i = 0; i < this->N; i++)
				out[i] = I * this->param_function_L[i];
		else if (this->model == "barotropic_triad_NL")
			for (std::size_t i = 0; i < this->N; i++)
				out[i] = 0.;
		else if (this->model == "SWE_triad")
			for (std::size_t i = 0; i < this->N; i++)
				out[i] = -I * this->param_function_L[i];
#endif
		else
			SWEETError("Unknown model " + this->model);

		return out;
	}

	ScalarDataArray function_L(
			ScalarDataArray &i_u,
			double i_dt,
			double i_sim_timestamp
		)
	{
		ScalarDataArray o_u = i_u;
		ScalarDataArray lL = this->lambda_L(i_dt, i_sim_timestamp);
		if (this->model == "ode1")
			for (std::size_t i = 0; i < o_u.number_of_elements; i++)
				o_u.set(i, std::sin(i_sim_timestamp));
		else
			o_u = lL * i_u;
			//for (std::size_t i = 0; i < o_u.number_of_elements; i++)
			//	o_u.set(i, lL.get(i) * i_u.get(i));
		return o_u;
	}

	ScalarDataArray function_N(
			ScalarDataArray &i_u,
			double i_dt,
			double i_sim_timestamp
		)
	{
		ScalarDataArray o_u = i_u;
		if (this->model == "ode1")
			for (std::size_t i = 0; i < o_u.number_of_elements; i++)
				o_u.set(i, std::sin(i_u.get(i)) * this->param_function_N[i]);
		else if (this->model == "cox_matthews_decay")
			for (std::size_t i = 0; i < o_u.number_of_elements; i++)
				o_u.set(i, std::sin(i_sim_timestamp) * this->param_function_N[i]);
#if SWEET_SCALAR_COMPLEX
		else if (this->model == "cox_matthews_oscillation")
			for (std::size_t i = 0; i < o_u.number_of_elements; i++)
				o_u.set(i, std::exp(I * i_sim_timestamp) * this->param_function_N[i]);
		else if (this->model == "barotropic_triad_NL")
		{
			double delta = param_function_extra[0] - param_function_extra[1] - param_function_extra[2];
			o_u.set(0, I * this->param_function_N[0] * i_u.get(1) * i_u.get(2)            * std::exp( I * delta * i_sim_timestamp));
			o_u.set(1, I * this->param_function_N[1] * i_u.get(0) * std::conj(i_u.get(2)) * std::exp(-I * delta * i_sim_timestamp));
			o_u.set(2, I * this->param_function_N[2] * i_u.get(0) * std::conj(i_u.get(1)) * std::exp(-I * delta * i_sim_timestamp));
		}
		else if (this->model == "SWE_triad")
		{
			o_u.set(0, I * this->param_function_N[0] * std::conj(i_u.get(1)) * i_u.get(2));
			o_u.set(1, I * this->param_function_N[1] * std::conj(i_u.get(0)) * i_u.get(2));
			o_u.set(2, I * this->param_function_N[2] *           i_u.get(0)  * i_u.get(1));
		}
#endif
		else
			SWEETError("Unknown model " + this->model);
		return o_u;
	}

public:
	virtual void run_timestep(
			///T &io_u,			///< prognostic variables
			ScalarDataArray &io_u,			///< prognostic variables

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
		ScalarDataArray u;
		u.setup(N_ode);
		for (int i = 0; i < N_ode; i++)
			u.set(i, io_data->get_pointer_to_data_Scalar()->simfields[i]);

		run_timestep(
				u,
				i_dt,
				i_sim_timestamp
			);

		for (int i = 0; i < N_ode; i++)
			io_data->get_pointer_to_data_Scalar()->simfields[i] = u.get(i);

	}

	// for parareal SL (not needed here)
	void set_previous_solution(
			Parareal_GenericData* i_data
	)
	{

		///u_prev = i_data->get_pointer_to_data_Scalar()->simfields[0];
		for (int i = 0; i < N_ode; i++)
			u_prev.set(i, i_data->get_pointer_to_data_Scalar()->simfields[i]);
	};
#endif

	void setup(
			std::string i_L,
			std::string i_N,
			std::string i_extra,
			std::string i_model
		)
	{

		this->u_prev.setup(N_ode);

		// get param L
		this->param_function_L = {};
		std::stringstream all_i_L = std::stringstream(i_L);
		while (all_i_L.good())
		{
			std::string str;
			getline(all_i_L, str, ',');
			this->param_function_L.push_back(atof(str.c_str()));
		}
		if ( this->param_function_L.size() != N_ode )
			SWEETError("param_function_L must contain N_ode values!");

		// get param N
		this->param_function_N = {};
		std::stringstream all_i_N = std::stringstream(i_N);
		while (all_i_N.good())
		{
			std::string str;
			getline(all_i_N, str, ',');
			this->param_function_N.push_back(atof(str.c_str()));
		}
		if ( this->param_function_N.size() != N_ode )
			SWEETError("param_function_N must contain N_ode values!");

		// get param extra
		this->param_function_extra = {};
		std::stringstream all_i_extra = std::stringstream(i_extra);
		while (all_i_extra.good())
		{
			std::string str;
			getline(all_i_extra, str, ',');
			this->param_function_extra.push_back(atof(str.c_str()));
		}

		this->model = i_model;

		this->N = N_ode;
	}

};


#endif
