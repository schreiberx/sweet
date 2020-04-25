/*
 * Adv_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_TIMESTEPPERS_HPP_

#include "Adv_Sphere_TS_interface.hpp"
#include "Adv_Sphere_TS_na_erk.hpp"
#include "Adv_Sphere_TS_na_sl.hpp"
#include "Adv_Sphere_TS_na_trajectories.hpp"



/**
 * SWE Sphere time steppers
 */
class Adv_Sphere_TimeSteppers
{
public:
	Adv_Sphere_TS_na_erk *na_erk = nullptr;
	Adv_Sphere_TS_na_sl *na_sl = nullptr;
	Adv_Sphere_TS_na_trajectories *na_trajectories = nullptr;
	Adv_Sphere_TS_interface *master = nullptr;



	Adv_Sphere_TimeSteppers()
	{
	}

	void reset()
	{
		if (na_erk != nullptr)
		{
			delete na_erk;
			na_erk = nullptr;
		}

		if (na_sl != nullptr)
		{
			delete na_sl;
			na_sl = nullptr;
		}

		if (na_trajectories != nullptr)
		{
			delete na_trajectories;
			na_trajectories = nullptr;
		}
	}


	void setup(
			const std::string &i_timestepping_method,
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	)
	{
		if (i_timestepping_method == "na_erk")
		{
			na_erk = new Adv_Sphere_TS_na_erk(i_simVars, i_op);
			na_erk->setup(i_simVars.disc.timestepping_order);

			master = &(Adv_Sphere_TS_interface&)*na_erk;
		}
		else if (i_timestepping_method == "na_sl")
		{
			na_sl = new Adv_Sphere_TS_na_sl(i_simVars, i_op);
			na_sl->setup(i_simVars.disc.timestepping_order);

			master = &(Adv_Sphere_TS_interface&)*na_sl;
		}
		else if (i_timestepping_method == "na_trajectories")
		{
			na_trajectories = new Adv_Sphere_TS_na_trajectories(i_simVars, i_op);
			na_trajectories->setup(i_simVars.disc.timestepping_order);

			master = &(Adv_Sphere_TS_interface&)*na_trajectories;
		}
		else
		{
			FatalError(std::string("No valid --timestepping-method provided: ")+i_timestepping_method);
		}
	}


	~Adv_Sphere_TimeSteppers()
	{
		reset();
	}
};




#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
