/*
 * ODE_Scalar_TimeSteppers.hpp
 *
 *  Created on: 08 Jun 2022
 *      Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_



#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_interface.hpp"

#include <sweet/SimulationVariables.hpp>

class ODE_Scalar_TimeSteppers
{
public:
	ODE_Scalar_TS_interface *master = nullptr;

	ODE_Scalar_TimeSteppers()
	{
	}

	void reset()
	{
		if (master != nullptr)
		{
			delete master;
			master = nullptr;
		}
	}

	void setup(
			//const std::string &i_timestepping_method,
			///int &i_timestepping_order,

			SimulationVariables &i_simVars
	)
	{
		reset();
		master = new ODE_Scalar_TS_interface;
		master->setup(
				atof(i_simVars.bogus.var[1].c_str()),
				atof(i_simVars.bogus.var[2].c_str())
			);
	}

	~ODE_Scalar_TimeSteppers()
	{
		reset();
	}
};


#endif
