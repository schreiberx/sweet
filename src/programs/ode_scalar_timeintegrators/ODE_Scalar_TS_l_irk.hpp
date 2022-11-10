/*
 * ODE_Scalar_TS_l_irk.hpp
 *
 *  Created on: 21 Oct 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_ODE_SCALAR_TS_L_IRK_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_ODE_SCALAR_TS_L_IRK_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_interface.hpp"


template <typename T>
class ODE_Scalar_TS_l_irk	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	int timestepping_order;

public:
	ODE_Scalar_TS_l_irk(
			SimulationVariables &i_simVars
		)
		: simVars(i_simVars)
	{
		///this->setup(simVars.disc.timestepping_order);
	}

	void setup(
			int i_order	///< order of RK time stepping method
	)
	{
		timestepping_order = i_order;

		if (timestepping_order != 2)
			SWEETError("ODE_Scalar_TS_l_irk: Only 2nd order IRK is supported. Please set --timestepping-order 2.");
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

		// Out vars
		ScalarDataArray u;

		double dt = i_dt;

		// Factor multiplying the implicit terms
		ScalarDataArray facI = 1. - i_dt * .5 * this->lambda_L(i_dt, i_simulation_timestamp);

		// Factor multiplying the explicit linear term
		ScalarDataArray facE = 1. + i_dt * .5 * this->lambda_L(i_dt, i_simulation_timestamp);

		// Solve
		u = 1. / facI * facE * io_u;

		// output data
		io_u = u;

	}

	virtual ~ODE_Scalar_TS_l_irk()
	{
	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
