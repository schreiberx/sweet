/*
 *  Created on: 07 Mar 2023
 *      Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#include "ODEScalarTS_ln_erk.hpp"


ODEScalarTS_ln_erk::ODEScalarTS_ln_erk()
{
}


ODEScalarTS_ln_erk::~ODEScalarTS_ln_erk()
{
}


/*
 * Setup
 */
bool ODEScalarTS_ln_erk::setup()
{
	this->param_a = shackODEScalarBenchmark->ode_parameters[0];
	this->param_b = shackODEScalarBenchmark->ode_parameters[1];
	return true;
}

void ODEScalarTS_ln_erk::runTimestep(
		double &io_u,		///< prognostic variables

		double i_dt,
		double i_simulation_timestamp
)
{
	assert(i_dt > 0);

	io_u += i_dt * (this->param_a * std::sin(io_u) + this->param_b * std::sin(i_simulation_timestamp));

}


