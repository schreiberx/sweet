/*
 * Burgers_Plane_TimeSteppers.hpp
 *
 *  Created on: 15 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TIMESTEPPERS_HPP_

#include "Burgers_Plane_TS_interface.hpp"

#include "Burgers_Plane_TS_l_irk_n_sl.hpp"
#include "Burgers_Plane_TS_l_irk_n_sl_forcing.hpp"
#include "Burgers_Plane_TS_ln_imex.hpp"
#include "Burgers_Plane_TS_ln_imex_forcing.hpp"
#include "Burgers_Plane_TS_ln_erk.hpp"
#include "Burgers_Plane_TS_ln_erk_forcing.hpp"
#include "Burgers_Plane_TS_ln_adomian.hpp"
#include "Burgers_Plane_TS_ln_cole_hopf.hpp"
#include "Burgers_Plane_TS_l_direct.hpp"
#include "Burgers_Plane_TS_l_erk.hpp"
#include "Burgers_Plane_TS_l_irk.hpp"



/**
 * Burgers Plane time steppers
 */
class Burgers_Plane_TimeSteppers
{
public:
	Burgers_Plane_TS_ln_erk *ln_erk = nullptr;
	Burgers_Plane_TS_ln_erk_forcing *ln_erk_forcing = nullptr;
	Burgers_Plane_TS_ln_imex *ln_imex = nullptr;
	Burgers_Plane_TS_ln_imex_forcing *ln_imex_forcing = nullptr;
	Burgers_Plane_TS_l_irk_n_sl *l_irk_n_sl = nullptr;
	Burgers_Plane_TS_l_irk_n_sl_forcing *l_irk_n_sl_forcing = nullptr;
	Burgers_Plane_TS_ln_adomian *ln_adomian = nullptr;
	Burgers_Plane_TS_ln_cole_hopf *ln_cole_hopf = nullptr;
	Burgers_Plane_TS_l_direct *l_direct = nullptr;
	Burgers_Plane_TS_l_erk *l_erk = nullptr;
	Burgers_Plane_TS_l_irk *l_irk = nullptr;

	Burgers_Plane_TS_interface *master = nullptr;

	Burgers_Plane_TimeSteppers()
	{
	}

	void reset()
	{
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		if (ln_erk_forcing != nullptr)
		{
			delete ln_erk_forcing;
			ln_erk_forcing = nullptr;
		}

		if (ln_imex != nullptr)
		{
			delete ln_imex;
			ln_imex = nullptr;
		}

		if (ln_imex_forcing != nullptr)
		{
			delete ln_imex_forcing;
			ln_imex_forcing = nullptr;
		}

		if (l_irk_n_sl != nullptr)
		{
			delete l_irk_n_sl;
			l_irk_n_sl = nullptr;
		}

		if (l_irk_n_sl_forcing != nullptr)
		{
			delete l_irk_n_sl_forcing;
			l_irk_n_sl_forcing = nullptr;
		}

		if (ln_adomian != nullptr)
		{
			delete ln_adomian;
			ln_adomian = nullptr;
		}

		if (ln_cole_hopf != nullptr)
		{
			delete ln_cole_hopf;
			ln_cole_hopf = nullptr;
		}

		if (l_direct != nullptr)
		{
			delete l_direct;
			l_direct = nullptr;
		}

		if (l_erk != nullptr)
		{
			delete l_erk;
			l_erk = nullptr;
		}

		if (l_irk != nullptr)
		{
			delete l_irk;
			l_irk = nullptr;
		}

	}



	void setup(
			const std::string &i_timestepping_method,
			int &i_timestepping_order,
			int &i_timestepping_order2,

			PlaneOperators &i_op,
			SimulationVariables &i_simVars
	)
	{
		reset();

		if (i_simVars.sim.CFL >= 0)
			FatalError("Only constant time step size supported with Burgers' equation, use negative CFL to set constant time step size");

		/// Always allocate analytical solution
		l_direct = new Burgers_Plane_TS_l_direct(i_simVars, i_op);
		l_direct->setup();
		ln_cole_hopf = new Burgers_Plane_TS_ln_cole_hopf(i_simVars, i_op);
		ln_cole_hopf->setup();

		if (i_timestepping_method == "ln_erk")
		{
			ln_erk = new Burgers_Plane_TS_ln_erk(i_simVars, i_op);
			ln_erk->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_erk;
		}
		else if (i_timestepping_method == "ln_erk_forcing")
		{
			ln_erk_forcing = new Burgers_Plane_TS_ln_erk_forcing(i_simVars, i_op);
			ln_erk_forcing->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_erk_forcing;
		}
		else if (i_timestepping_method == "ln_imex")
		{
			ln_imex= new Burgers_Plane_TS_ln_imex(i_simVars, i_op);
			ln_imex->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_imex;
		}
		else if (i_timestepping_method == "ln_imex_forcing")
		{
			ln_imex_forcing= new Burgers_Plane_TS_ln_imex_forcing(i_simVars, i_op);
			ln_imex_forcing->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_imex_forcing;
		}
		else if (i_timestepping_method == "l_irk_n_sl")
		{
			l_irk_n_sl = new Burgers_Plane_TS_l_irk_n_sl(i_simVars, i_op);
			l_irk_n_sl->setup();

			master = &(Burgers_Plane_TS_interface&)*l_irk_n_sl;
		}
		else if (i_timestepping_method == "l_irk_n_sl_forcing")
		{
			l_irk_n_sl_forcing = new Burgers_Plane_TS_l_irk_n_sl_forcing(i_simVars, i_op);
			l_irk_n_sl_forcing->setup();

			master = &(Burgers_Plane_TS_interface&)*l_irk_n_sl_forcing;
		}
		else if (i_timestepping_method == "ln_adomian")
		{
			ln_adomian= new Burgers_Plane_TS_ln_adomian(i_simVars, i_op);
			ln_adomian->setup();

			master = &(Burgers_Plane_TS_interface&)*ln_adomian;
		}
		else if (i_timestepping_method == "ln_cole_hopf")
		{
			ln_cole_hopf= new Burgers_Plane_TS_ln_cole_hopf(i_simVars, i_op);
			ln_cole_hopf->setup();

			master = &(Burgers_Plane_TS_interface&)*ln_cole_hopf;
		}
		else if (i_timestepping_method == "l_direct")
		{
			l_direct= new Burgers_Plane_TS_l_direct(i_simVars, i_op);
			l_direct->setup();

			master = &(Burgers_Plane_TS_interface&)*l_direct;
		}
		else if (i_timestepping_method == "l_erk")
		{
			l_erk= new Burgers_Plane_TS_l_erk(i_simVars, i_op);
			l_erk->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*l_erk;
		}
		else if (i_timestepping_method == "l_irk")
		{
			l_irk= new Burgers_Plane_TS_l_irk(i_simVars, i_op);
			l_irk->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*l_irk;
		}
		//
		else
		{
			std::cout << "Unknown method: " << i_timestepping_method << std::endl;
			std::cout << "Available --timestepping-method :"  << std::endl;
			std::cout << "      l_direct           : Linear: analytical solution to diffusion operator"  << std::endl;
			std::cout << "      l_erk              : Linear: explicit RK scheme"  << std::endl;
			std::cout << "      l_irk              : Linear: implicit RK scheme"  << std::endl;
			std::cout << "      ln_cole_hopf       : Non-linear: analytic solution to Burgers' equation"  << std::endl;
			std::cout << "      ln_erk             : Non-linear: explicit RK scheme"  << std::endl;
			std::cout << "      ln_erk_forcing     : Non-linear: explicit RK scheme with forcing term"  << std::endl;
			std::cout << "      ln_imex            : Non-linear: implicit-explicit RK scheme"  << std::endl;
			std::cout << "      ln_imex_forcing    : Non-linear: implicit-explicit RK scheme with forcing term"  << std::endl;
			std::cout << "      l_irk_n_sl         : Non-linear: implicit RK on semi-Lagrangian formulation"  << std::endl;
			std::cout << "      l_irk_n_sl_forcing : Non-linear: implicit RK on semi-Lagrangian formulation with forcing"  << std::endl;
			std::cout << "      ln_adomian         : Non-linear: Adomian decomposition method"  << std::endl;
			FatalError("No valid --timestepping-method provided");
		}
	}



	~Burgers_Plane_TimeSteppers()
	{
		reset();
	}
};




#endif /* SRC_PROGRAMS_BURGERS_PLANE_TIMESTEPPERS_HPP_ */
