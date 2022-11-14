/*
 * ODE_Scalar_TS_l_cn_n_settls.hpp
 *
 *  Created on: 04 Oct 2022
 *      Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 *  Modif on: 11 July 2017
 *  	Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Name: l_cn_na_sl_nd_settls
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TS_L_CN_N_SETTLS_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TS_L_CN_N_SETTLS_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_interface.hpp"


template <typename T>
class ODE_Scalar_TS_l_cn_n_settls	: public ODE_Scalar_TS_interface<T>
{
	SimulationVariables &simVars;

	////ScalarDataArray u_prev;

public:
	ODE_Scalar_TS_l_cn_n_settls(
			SimulationVariables &i_simVars
		)
		:
		simVars(i_simVars)
	{
	}

	void setup(
	)
	{
	}

	void run_timestep(
			///T &io_u,	///< prognostic variables
			ScalarDataArray &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	)
	{

		if (i_dt <= 0)
			SWEETError("SWE_Plane_TS_l_cn_na_sl_nd_settls: Only constant time step size allowed (Please set --dt)");

		if (i_simulation_timestamp == 0)
		{
#if (!SWEET_PARAREAL) && (!SWEET_XBRAID)
			/*
			 * First time step
			 */
			this->u_prev = io_u;
#endif
		}


		if (std::abs(i_dt - 0.08)< 1e-10 && i_simulation_timestamp < 0.25)
		{
		std::cout << "TIME " << i_simulation_timestamp << " " << i_dt << std::endl;
		std::cout << "u_prev " << this->u_prev[0] << " " << this->u_prev[1] << " " << this->u_prev[2] << std::endl;
		std::cout << "u " << io_u[0] << " " << io_u[1] << " " << io_u[2] << std::endl;
		}

		// Out vars
		ScalarDataArray u;
		u.setup(N_ode);

		double dt = i_dt;

		// Nonlinear term at t_n
		ScalarDataArray nonlin = this->function_N(io_u, i_dt, i_simulation_timestamp);

		// Nonlinear term at t_{n-1}
		ScalarDataArray nonlin1 = this->function_N(this->u_prev, i_dt, i_simulation_timestamp);

		// Extrapolate + average nonlinear term
		nonlin = i_dt * .5 * (2. * nonlin - nonlin1 + nonlin);

		// Factor multiplying the implicit terms
		ScalarDataArray facI = 1. - i_dt * .5 * this->lambda_L(i_dt, i_simulation_timestamp);

		// Factor multiplying the explicit linear term
		ScalarDataArray facE = 1. + i_dt * .5 * this->lambda_L(i_dt, i_simulation_timestamp);

		// Solve
		u = 1. / facI * (facE * io_u + nonlin);

		this->u_prev = io_u;

		// output data
		io_u = u;

		if (std::abs(i_dt - 0.08)< 1e-5 && i_simulation_timestamp < 0.25)
		{
		std::cout << "new u " << io_u[0] << " " << io_u[1] << " " << io_u[2] << std::endl;
		}

	}

	virtual ~ODE_Scalar_TS_l_cn_n_settls()
	{
	}
};

#endif /* SRC_PROGRAMS_ODE_SCALAR_TS_L_CN_N_SETTLS_HPP_ */
