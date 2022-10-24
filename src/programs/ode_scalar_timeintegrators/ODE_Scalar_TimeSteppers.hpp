/*
 * ODE_Scalar_TimeSteppers.hpp
 *
 *  Created on: 08 Jun 2022
 *      Author: Joao Steinstraessrt <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_

#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_interface.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_ln_erk.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_l_irk.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_n_erk.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_l_irk_n_erk.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_l_cn_n_settls.hpp"
#include "../ode_scalar_timeintegrators/ODE_Scalar_TS_l_exp_n_etdrk.hpp"

#include <sweet/SimulationVariables.hpp>

template <typename T>
class ODE_Scalar_TimeSteppers
{
public:
	ODE_Scalar_TS_ln_erk<T> *ln_erk = nullptr;
	ODE_Scalar_TS_n_erk<T> *n_erk = nullptr;
	ODE_Scalar_TS_l_irk<T> *l_irk = nullptr;
	ODE_Scalar_TS_l_irk_n_erk<T> *l_irk_n_erk = nullptr;
	ODE_Scalar_TS_l_cn_n_settls<T> *l_cn_n_settls = nullptr;
	ODE_Scalar_TS_l_exp_n_etdrk<T> *l_exp_n_etdrk = nullptr;

	ODE_Scalar_TS_interface<T> *master = nullptr;

	ODE_Scalar_TimeSteppers()
	{
	}

	void reset()
	{
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		if (n_erk != nullptr)
		{
			delete n_erk;
			n_erk = nullptr;
		}

		if (l_irk != nullptr)
		{
			delete l_irk;
			l_irk = nullptr;
		}

		if (l_irk_n_erk != nullptr)
		{
			delete l_irk_n_erk;
			l_irk_n_erk = nullptr;
		}

		if (l_cn_n_settls != nullptr)
		{
			delete l_cn_n_settls;
			l_cn_n_settls = nullptr;
		}

		if (l_exp_n_etdrk != nullptr)
		{
			delete l_exp_n_etdrk;
			l_exp_n_etdrk = nullptr;
		}

		///////if (master != nullptr)
		///////{
		///////	delete master;
		///////	master = nullptr;
		///////}
	}

	void setup(
			const std::string &i_timestepping_method,
			int &i_timestepping_order,
			SimulationVariables &i_simVars
	)
	{
		reset();

		if (i_timestepping_method == "ln_erk")
		{
			ln_erk = new ODE_Scalar_TS_ln_erk<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*ln_erk;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			ln_erk->setup(
						i_timestepping_order
					);

		}

		else if (i_timestepping_method == "n_erk")
		{
			n_erk = new ODE_Scalar_TS_n_erk<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*n_erk;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			n_erk->setup(
						i_timestepping_order
					);

		}

		else if (i_timestepping_method == "l_irk")
		{
			l_irk = new ODE_Scalar_TS_l_irk<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*l_irk;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			l_irk->setup(
						i_timestepping_order
					);

		}

		else if (i_timestepping_method == "l_irk_n_erk")
		{
			l_irk_n_erk = new ODE_Scalar_TS_l_irk_n_erk<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*l_irk_n_erk;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			l_irk_n_erk->setup(
						i_timestepping_order
					);

		}


		else if (i_timestepping_method == "l_cn_n_settls")
		{
			l_cn_n_settls = new ODE_Scalar_TS_l_cn_n_settls<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*l_cn_n_settls;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			l_cn_n_settls->setup(
					);

		}

		else if (i_timestepping_method == "l_exp_n_etdrk")
		{
			l_exp_n_etdrk = new ODE_Scalar_TS_l_exp_n_etdrk<T>(i_simVars);
			master = &(ODE_Scalar_TS_interface<T>&)*l_exp_n_etdrk;

			master->setup(
					i_simVars.bogus.var[3],
					i_simVars.bogus.var[4],
					i_simVars.bogus.var[5],
					i_simVars.bogus.var[6]
				);

			l_exp_n_etdrk->setup(
						i_simVars.rexi,
						i_timestepping_order
					);

		}

		else //Help menu with list of schemes
		{
			std::cout << "Unknown method: " << i_timestepping_method << std::endl;
			std::cout << "Available --timestepping-method :"  << std::endl;
			std::cout << "      l_direct       : Linear:     Analytical solution to linear SW operator"  << std::endl;
			std::cout << "      l_erk          : Linear:     Explicit RK scheme (supports FD-C staggering)"  << std::endl;
			std::cout << "      l_erk_n_erk    : Non-linear: Linear RK, Non-linear RK, Strang-split"  << std::endl;
			std::cout << "      l_cn           : Linear:     Crank-Nicolson (CN) scheme"  << std::endl;
			std::cout << "      l_cn_n_erk     : Non-linear: Linear CN, Non-linear RK, Strang-split"<< std::endl;
			std::cout << "      l_rexi         : Linear:     Pure REXI, our little dog."<< std::endl;
			std::cout << "      l_rexi_n_erk   : Non-linear: Linear REXI, Non-linear RK, Strang-split"<< std::endl;
			std::cout << "      l_irk          : Linear:     Implicit Euler"  << std::endl;
			std::cout << "      l_irk_n_erk    : Non-linear: Linear Implicit Euler, Non-linear RK, Strang-split"  << std::endl;
			std::cout << "      ln_erk         : Non-linear: Linear and nonlinear solved jointly with erk (supports FD-C staggering)"  << std::endl;
			std::cout << "      l_rexi_na_sl_nd_settls   : Non-linear: Linear Rexi, Advection: Semi-Lag, Nonlinear-diverg: SETTLS"  << std::endl;
			std::cout << "      l_cn_na_sl_nd_settls     : Non-linear: Linear CN, Advection: Semi-Lag, Nonlinear-diverg: SETTLS"  << std::endl;
			std::cout << "      l_rexi_n_etdrk           : Non-linear: Linear REXI, Non-linear: ETDRK"  << std::endl;
			std::cout << "      l_rexi_na_sl_nd_etdrk    : Non-linear: Linear REXI, Advection: Semi-Lag, Nonlinear-diverg: ETDRK"  << std::endl;

			SWEETError("No valid --timestepping-method provided");
		}




		////////master = new ODE_Scalar_TS_interface;
		////////master->setup(
		////////		atof(i_simVars.bogus.var[1].c_str()),
		////////		atof(i_simVars.bogus.var[2].c_str())
		////////	);
	}

	~ODE_Scalar_TimeSteppers()
	{
		reset();
	}
};


#endif
