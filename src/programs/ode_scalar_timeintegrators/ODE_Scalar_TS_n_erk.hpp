/*
 * ODE_Scalar_TS_n_erk.hpp
 *
 *  Created on: 30 Sep 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_N_ERK_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_N_ERK_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "ODE_Scalar_TS_interface.hpp"


template <typename T>
class ODE_Scalar_TS_n_erk	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	int timestepping_order;

public:
	ODE_Scalar_TS_n_erk(
			SimulationVariables &i_simVars
		)
		:
		simVars(i_simVars)
	{
		setup(simVars.disc.timestepping_order);
	}

	void setup(
			int i_order	///< order of RK time stepping method
	)
	{
		timestepping_order = i_order;
	}

	void setup(
			std::string i_L,
			std::string i_N,
			std::string i_extra,
			std::string i_model
		)
	{
		ODE_Scalar_TS_interface<T>::setup(i_L, i_N, i_extra, i_model);
	}

	void run_timestep(
			ScalarDataArray &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{
		if (i_dt <= 0)
			SWEETError("ODE_Scalar_TS_ln_erk: Only constant time step size allowed");

		if (this->timestepping_order == 1)
			io_u += i_dt * (this->function_N(io_u, i_dt, i_simulation_timestamp));
		else if (this->timestepping_order == 2)
		{
			ScalarDataArray k1 = this->function_N(io_u, i_dt, i_simulation_timestamp);
			ScalarDataArray u2 = io_u + i_dt * k1;
			ScalarDataArray k2 = this->function_N(u2, i_dt, i_simulation_timestamp + i_dt);
			io_u += .5 * i_dt * (k1 + k2);
		}
	}



	virtual ~ODE_Scalar_TS_n_erk()
	{
	}
};

#endif /* SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_N_ERK_HPP_ */
