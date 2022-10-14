/*
 * ODE_Scalar_TS_l_exp_n_etdrk.hpp
 *
 *  Created on: 28 September 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_EXP_N_ETDRK_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_EXP_N_ETDRK_HPP_

#include <limits>
#include "ODE_Scalar_TS_interface.hpp"
#include "ODE_Scalar_TS_l_exp.hpp"

template <typename T>
class ODE_Scalar_TS_l_exp_n_etdrk	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	ODE_Scalar_TS_l_exp<T> ts_phi0_exp;
	ODE_Scalar_TS_l_exp<T> ts_phi1_exp;
	ODE_Scalar_TS_l_exp<T> ts_phi2_exp;

	ODE_Scalar_TS_l_exp<T> ts_ups0_exp;
	ODE_Scalar_TS_l_exp<T> ts_ups1_exp;
	ODE_Scalar_TS_l_exp<T> ts_ups2_exp;
	ODE_Scalar_TS_l_exp<T> ts_ups3_exp;

	int timestepping_order;


public:
	ODE_Scalar_TS_l_exp_n_etdrk(
			SimulationVariables &i_simVars
		)
		:
		simVars(i_simVars),
		ts_phi0_exp(simVars),
		ts_phi1_exp(simVars),
		ts_phi2_exp(simVars),

		ts_ups0_exp(simVars),
		ts_ups1_exp(simVars),
		ts_ups2_exp(simVars),
		ts_ups3_exp(simVars)
	{
	}

	void setup(
			EXP_SimulationVariables &i_expSimVars,
			int i_timestepping_order
	)
	{
		timestepping_order = i_timestepping_order;

		if (timestepping_order == 1)
		{
			ts_phi0_exp.setup(i_expSimVars, "phi0", simVars.timecontrol.current_timestep_size);
			ts_phi1_exp.setup(i_expSimVars, "phi1", simVars.timecontrol.current_timestep_size);

			ts_phi0_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_phi1_exp.setup(this->param_function_L, this->param_function_N, this->model);
		}
		else if (timestepping_order == 2)
		{
			ts_phi0_exp.setup(i_expSimVars, "phi0", simVars.timecontrol.current_timestep_size);
			ts_phi1_exp.setup(i_expSimVars, "phi1", simVars.timecontrol.current_timestep_size);
			ts_phi2_exp.setup(i_expSimVars, "phi2", simVars.timecontrol.current_timestep_size);

			ts_phi0_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_phi1_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_phi2_exp.setup(this->param_function_L, this->param_function_N, this->model);
		}
		else if (timestepping_order == 4)
		{
			ts_phi0_exp.setup(i_expSimVars, "phi0", simVars.timecontrol.current_timestep_size*0.5);
			ts_phi1_exp.setup(i_expSimVars, "phi1", simVars.timecontrol.current_timestep_size*0.5);
			ts_phi2_exp.setup(i_expSimVars, "phi2", simVars.timecontrol.current_timestep_size*0.5);

			ts_ups0_exp.setup(i_expSimVars, "phi0", simVars.timecontrol.current_timestep_size);
			ts_ups1_exp.setup(i_expSimVars, "ups1", simVars.timecontrol.current_timestep_size);
			ts_ups2_exp.setup(i_expSimVars, "ups2", simVars.timecontrol.current_timestep_size);
			ts_ups3_exp.setup(i_expSimVars, "ups3", simVars.timecontrol.current_timestep_size);

			ts_phi0_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_phi1_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_phi2_exp.setup(this->param_function_L, this->param_function_N, this->model);

			ts_ups0_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_ups1_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_ups2_exp.setup(this->param_function_L, this->param_function_N, this->model);
			ts_ups3_exp.setup(this->param_function_L, this->param_function_N, this->model);
		}
	}

	void euler_timestep_update_nonlinear(
			T &i_u,	///< prognostic variables

			T &o_u,	///< time updates

			double i_dt,
			double i_simulation_timestamp
	)
	{
		if (i_dt <= 0)
			SWEETError("ODE_Scalar_TS_ln_erk: Only constant time step size allowed");

		o_u = this->function_N(i_u, i_dt, i_simulation_timestamp);

	}

	void run_timestep(
			T &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{
		if (i_dt <= 0)
			SWEETError("ODE_Scalar_TS_l_phi0_n_edt: Only constant time step size allowed");

		if (timestepping_order == 1)
		{
			/*
			 * U_{1} = \psi_{0}( \Delta t L ) U_{0}
			 * 			+\Delta t \psi_{1}(\Delta tL) N(U_{0}).
			 */

			T phi0_Un_u;
			ts_phi0_exp.run_timestep(
					io_u,
					phi0_Un_u,
					i_dt,
					i_simulation_timestamp
				);

			T FUn_u;
			euler_timestep_update_nonlinear(
					io_u,
					FUn_u,
					i_dt,
					i_simulation_timestamp
			);

			T phi1_FUn_u;

			ts_phi1_exp.run_timestep(
					FUn_u,
					phi1_FUn_u,
					i_dt,
					i_simulation_timestamp
				);

			io_u = phi0_Un_u + i_dt*phi1_FUn_u;
		}
		else if (timestepping_order == 2)
		{

			/*
			 * A_{n}=\psi_{0}(\Delta tL)U_{n}+\Delta t\psi_{1}(\Delta tL)F(U_{n})
			 */

			T phi0_Un_u;

			ts_phi0_exp.run_timestep(
					io_u,
					phi0_Un_u,
					i_dt,
					i_simulation_timestamp
				);

			T FUn_u;

			euler_timestep_update_nonlinear(
					io_u,
					FUn_u,
					i_dt,
					i_simulation_timestamp
			);

			T phi1_FUn_u;

			ts_phi1_exp.run_timestep(
					FUn_u,
					phi1_FUn_u,
					i_dt,
					i_simulation_timestamp
				);

			T A_u = phi0_Un_u + i_dt*phi1_FUn_u;

			/*
			 * U_{n+1} = A_{n}+ \Delta t \psi_{2}(\Delta tL)
			 * 				\left(F(A_{n},t_{n}+\Delta t)-F(U_{n})\right)
			 */

			T FAn_u;

			euler_timestep_update_nonlinear(
					A_u,
					FAn_u,
					i_dt,
					i_simulation_timestamp
			);

			T phi2_X_u;

			ts_phi2_exp.run_timestep(
					FAn_u - FUn_u,

					phi2_X_u,

					i_dt,
					i_simulation_timestamp
				);

			io_u = A_u + i_dt*phi2_X_u;
		}
		else if (timestepping_order == 4)
		{
			double dt = i_dt;
			double dt_half = dt*0.5;

			/*
			 * Precompute commonly used terms
			 */
			T phi0_Un_u;

			ts_phi0_exp.run_timestep(
					io_u,
					phi0_Un_u,
					dt_half,
					i_simulation_timestamp
				);

			T FUn_u;

			euler_timestep_update_nonlinear(
					io_u,
					FUn_u,
					i_dt,
					i_simulation_timestamp
			);

			/*
			 * Some commonly shared buffers
			 */

			T phi1_u;

			/*
			 * A_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + \Delta t\psi_{1}(0.5*\Delta tL) F(U_{n})
			 */
			ts_phi1_exp.run_timestep(
					FUn_u,
					phi1_u,
					dt_half,
					i_simulation_timestamp
				);

			T A_u = phi0_Un_u + dt_half*phi1_u;


			/*
			 * B_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5*\Delta tL) F(A_{n}, t_{n} + 0.5*\Delta t)
			 */

			T FAn_u;

			euler_timestep_update_nonlinear(
					A_u,
					FAn_u,
					i_dt,
					i_simulation_timestamp + dt_half
			);

			ts_phi1_exp.run_timestep(
					FAn_u,
					phi1_u,
					dt_half,
					i_simulation_timestamp
				);

			T B_u = phi0_Un_u + dt_half*phi1_u;


			/*
			 * C_{n} = \psi_{0}(0.5*\Delta tL)U_{n} + 0.5*\Delta t\psi_{1}(0.5* \Delta tL) ( 2 F(B_{n},t_{n} + 0.5*\Delta t)-F(U_{n},t_{n})).
			 */

			T phi0_An_u;

			ts_phi0_exp.run_timestep(
					A_u,
					phi0_An_u,
					dt_half,
					i_simulation_timestamp
				);

			T FBn_u;

			euler_timestep_update_nonlinear(
					B_u,
					FBn_u,
					i_dt,
					i_simulation_timestamp + dt_half
			);

			ts_phi1_exp.run_timestep(
					2.0*FBn_u - FUn_u,
					phi1_u,
					dt_half,
					i_simulation_timestamp
				);

			T C_u = phi0_An_u + dt_half*phi1_u;


			/*
			 * R0 - R3
			 */
			T FCn_u;

			euler_timestep_update_nonlinear(
					C_u,
					FCn_u,
					i_dt,
					i_simulation_timestamp + dt
			);

			T R0_u = io_u;

			T &R1_u = FUn_u;

			T R2_u = FAn_u + FBn_u;

			T &R3_u = FCn_u;


			/*
			 * U_{n+1} =
			 * 		\psi_{0}(\Delta tL)R_{0}
			 * 			+ \Delta t
			 * 			(
			 * 				  \upsilon_{1}(\Delta tL) R_{1} +
			 * 				2*\upsilon_{2}(\Delta tL) R_{2} +
			 * 				  \upsilon_{3}(\Delta tL) R_{3}
			 * 			)
			 */
			ts_phi0_exp.run_timestep(
					R0_u,
					dt,
					i_simulation_timestamp
				);

			ts_ups1_exp.run_timestep(
					R1_u,
					dt,
					i_simulation_timestamp
				);

			ts_ups2_exp.run_timestep(
					R2_u,
					dt,
					i_simulation_timestamp
				);

			ts_ups3_exp.run_timestep(
					R3_u,
					dt,
					i_simulation_timestamp
				);

			io_u = R0_u + dt*(R1_u + 2.0*R2_u + R3_u);
		}
		else
		{
			SWEETError("TODO: This order is not implemented, yet!");
		}

	}



	virtual ~ODE_Scalar_TS_l_exp_n_etdrk()
	{
	}
};

#endif /* SRC_PROGRAMS_ODE_SCALAR_TIMEINTEGRATORS_ODE_SCALAR_TS_L_EXP_N_ETDRK_HPP_ */
