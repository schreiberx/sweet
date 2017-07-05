/*
 * Burgers_Plane_TimeSteppers.hpp
 *
 *  Created on: 15 June 2017
 *      Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_Burgers_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_Burgers_PLANE_TIMESTEPPERS_HPP_

#include "Burgers_Plane_TS_interface.hpp"

#include "Burgers_Plane_TS_l_irk_n_sl.hpp"
#include "Burgers_Plane_TS_ln_imex.hpp"
#include "Burgers_Plane_TS_ln_erk.hpp"



/**
 * Burgers Plane time steppers
 */
class Burgers_Plane_TimeSteppers
{
public:
	Burgers_Plane_TS_ln_erk *ln_erk = nullptr;
	Burgers_Plane_TS_ln_imex *ln_imex = nullptr;
	Burgers_Plane_TS_l_irk_n_sl *l_irk_n_sl = nullptr;

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

		if (ln_imex != nullptr)
		{
			delete ln_imex;
			ln_imex = nullptr;
		}

		if (l_irk_n_sl != nullptr)
		{
			delete l_irk_n_sl;
			l_irk_n_sl = nullptr;
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

		/*
		/// Always allocate analytical solution
		l_direct = new Burgers_Plane_TS_l_direct(i_simVars, i_op);
		*/

		if (i_timestepping_method == "ln_erk")
		{
			ln_erk = new Burgers_Plane_TS_ln_erk(i_simVars, i_op);
			ln_erk->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_erk;
		}
		else if (i_timestepping_method == "ln_imex")
		{
			ln_imex= new Burgers_Plane_TS_ln_imex(i_simVars, i_op);
			ln_imex->setup(i_timestepping_order);

			master = &(Burgers_Plane_TS_interface&)*ln_imex;
		}
		else if (i_timestepping_method == "l_irk_n_sl")
		{
			l_irk_n_sl = new Burgers_Plane_TS_l_irk_n_sl(i_simVars, i_op);
			l_irk_n_sl->setup();

			master = &(Burgers_Plane_TS_interface&)*l_irk_n_sl;
		}
		//
		else
		{
			std::cout << i_timestepping_method << std::endl;
			FatalError("No valid --timestepping-method provided");
		}
	}



	~Burgers_Plane_TimeSteppers()
	{
		reset();
	}
};




#endif /* SRC_PROGRAMS_Burgers_PLANE_TIMESTEPPERS_HPP_ */
