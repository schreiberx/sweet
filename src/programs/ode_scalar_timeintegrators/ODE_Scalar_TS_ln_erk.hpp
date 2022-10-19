/*
 * ODE_Scalar_TS_ln_erk.hpp
 *
 *  Created on: 30 Sep 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_LN_ERK_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_LN_ERK_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "ODE_Scalar_TS_interface.hpp"


template <typename T>
class ODE_Scalar_TS_ln_erk	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	int timestepping_order;

public:
	ODE_Scalar_TS_ln_erk(
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

	void run_timestep(
			///T &io_u,	///< prognostic variables
			ScalarDataArray &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{
		if (i_dt <= 0)
			SWEETError("ODE_Scalar_TS_ln_erk: Only constant time step size allowed");

		/////io_u += i_dt * (this->param_function_L * std::sin(io_u) + this->param_function_N * std::sin(i_simulation_timestamp));
		io_u += i_dt * (  this->function_L(io_u, i_dt, i_simulation_timestamp)
				+ this->function_N(io_u, i_dt, i_simulation_timestamp));
	}



	virtual ~ODE_Scalar_TS_ln_erk()
	{
	}
};

#endif /* SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_LN_ERK_HPP_ */
