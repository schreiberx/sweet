/*
 * Burgers_Plane_TS_ln_imex_init_sl_forcing.cpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_ln_imex_init_sl_forcing.hpp"


void Burgers_Plane_TS_ln_imex_init_sl_forcing::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_ln_imex_init_sl_forcing: Only constant time step size allowed");


	double dt = i_fixed_dt;

	Staggering staggering;
	assert(staggering.staggering_type == 'a');

	if (simVars.bogus.var[0] == 0)
	{
		std::cout << "Bin in SL" << std::endl;
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk_n_sl_forcing.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				dt,
				i_simulation_timestamp
		);
	}
	else
	{
		std::cout << "Bin in IMEX" << std::endl;
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_ln_imex_forcing.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				dt,
				i_simulation_timestamp
		);
	}

}



/*
 * Setup
 */
void Burgers_Plane_TS_ln_imex_init_sl_forcing::setup()
{
	ts_l_irk_n_sl_forcing.setup();
	ts_ln_imex_forcing.setup(simVars.disc.timestepping_order);

}


Burgers_Plane_TS_ln_imex_init_sl_forcing::Burgers_Plane_TS_ln_imex_init_sl_forcing(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		ts_l_irk_n_sl_forcing(simVars, op),
		ts_ln_imex_forcing(simVars, op)
{
}



Burgers_Plane_TS_ln_imex_init_sl_forcing::~Burgers_Plane_TS_ln_imex_init_sl_forcing()
{
}

