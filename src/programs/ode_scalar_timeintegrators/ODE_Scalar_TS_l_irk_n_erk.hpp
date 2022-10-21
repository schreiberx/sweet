/*
 * ODE_Scalar_TS_l_irk_n_erk.hpp
 *
 *  Created on: 21 Oct 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_IRK_N_ERK_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_IRK_N_ERK_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include "ODE_Scalar_TS_interface.hpp"
#include "ODE_Scalar_TS_l_irk.hpp"


template <typename T>
class ODE_Scalar_TS_l_irk_n_erk	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	int timestepping_order;

	ODE_Scalar_TS_l_irk<T> ts_l_irk;
	ODE_Scalar_TS_n_erk<T> ts_n_erk;


public:
	ODE_Scalar_TS_l_irk_n_erk(
			SimulationVariables &i_simVars
		):
		simVars(i_simVars),
		ts_l_irk(simVars),
		ts_n_erk(simVars)
	{
		setup(simVars.disc.timestepping_order);
	}

	void setup(
			int i_order
	)
	{
		this->timestepping_order = i_order;

		ts_l_irk.setup(2); // CHECK THIS
		if (this->timestepping_order == 1) 
			ts_n_erk.setup(1); // CHECK THIS
		else if (this->timestepping_order == 2)
			ts_n_erk.setup(2); // CHECK THIS
		else
			SWEETError("Unsupported order.");

		////ODE_Scalar_TS_interface<T> *master = nullptr;
		////master = &(ODE_Scalar_TS_interface<T>&)*ts_l_irk;
		////master->setup(this->param_function_L, this->param_function_N, this->param_function_extra, this->model);
		ts_l_irk.setup(
					simVars.bogus.var[3],
					simVars.bogus.var[4],
					simVars.bogus.var[5],
					simVars.bogus.var[6]
		)
	}

	void run_timestep(
			ScalarDataArray &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{
		if (timestepping_order == 1)
		{
			// first order IRK for linear
			ts_l_irk.run_timestep(
					io_u,
					i_dt,
					i_simulation_timestamp
				);

			// first order explicit for non-linear
			ts_n_erk.run_timestep(
					io_u,
					i_dt,
					i_simulation_timestamp
				);
		}
		else if (timestepping_order == 2)
		{
			// HALF time step for linear part
			ts_l_irk.run_timestep(
					io_u,
					i_dt*0.5,
					i_simulation_timestamp
				);

			// FULL time step for non-linear part
			ts_n_erk.run_timestep(
					io_u,
					i_dt,
					i_simulation_timestamp
				);

			// HALF time step for linear part
			ts_l_irk.run_timestep(
					io_u,
					i_dt*0.5,
					i_simulation_timestamp+i_dt*0.5	/* TODO: CHECK THIS, THIS MIGHT BE WRONG!!! */
				);
		}
		else
		{
			SWEETError("Not yet supported!");
		}

	}

	ODE_Scalar_TS_l_irk<T>& get_implicit_timestepper()
	{
		return ts_l_irk;
	}

	virtual ~ODE_Scalar_TS_l_irk_n_erk()
	{
	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
