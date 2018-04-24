/*
 * Burgers_Plane_TS_l_irk_n_adomian.cpp
 *
 *  Created on: 14 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 *
 */

#include "Burgers_Plane_TS_l_irk_n_adomian.hpp"


void Burgers_Plane_TS_l_irk_n_adomian::run_timestep(
		PlaneData &io_u,	///< prognostic variables
		PlaneData &io_v,	///< prognostic variables
		PlaneData &io_u_prev,	///< prognostic variables
		PlaneData &io_v_prev,	///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp
)
{
	if (i_fixed_dt <= 0)
		FatalError("Burgers_Plane_TS_l_irk_n_adomian: Only constant time step size allowed");

	double dt = i_fixed_dt;

	if (simVars.disc.timestepping_order == 2)
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				0.5*dt,
				i_simulation_timestamp
		);
	}

	ts_n_adomian.run_timestep(
			io_u, io_v,
			io_u_prev, io_v_prev,
			dt,
			i_simulation_timestamp
	);

	if (simVars.disc.timestepping_order == 2)
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk.run_timestep(
				io_u, io_v,
				io_u_prev, io_v_prev,
				0.5*dt,
				i_simulation_timestamp
		);
	}
	else
	{
		// Run implicit Runge-Kutta on Burgers' equation in SL form
		ts_l_irk.run_timestep(
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
void Burgers_Plane_TS_l_irk_n_adomian::setup()
{
	ts_l_irk.setup(simVars.disc.timestepping_order);
	ts_n_adomian.setup(simVars.disc.timestepping_order2);
}


Burgers_Plane_TS_l_irk_n_adomian::Burgers_Plane_TS_l_irk_n_adomian(
		SimulationVariables &i_simVars,
		PlaneOperators &i_op
)	:
		simVars(i_simVars),
		op(i_op),

		ts_l_irk(simVars, op),
		ts_n_adomian(simVars,op)
{
}



Burgers_Plane_TS_l_irk_n_adomian::~Burgers_Plane_TS_l_irk_n_adomian()
{
}

