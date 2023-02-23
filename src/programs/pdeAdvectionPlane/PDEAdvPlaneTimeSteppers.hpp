/*
 * Adv_Plane_TimeSteppers.hpp
 *
 *  Created on: 4th April 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTION_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_PDE_ADVECTION_PLANE_TIMESTEPPERS_HPP_

#include "PDEAdvPlaneTS_interface.hpp"
#include "PDEAdvPlaneTS_na_erk.hpp"
//#include "PDEAdvPlaneTS_na_sl.hpp"

#include <sweet/ErrorBase.hpp>
#include <sweet/shacks/ShackDictionary.hpp>
#include "ShackPDEAdvectionPlaneTimeDiscretization.hpp"

/**
 * Advection for the plane time integrators
 */
class PDEAdvPlaneTimeSteppers
{
public:
	sweet::ErrorBase error;

	PDEAdvPlaneTS_na_erk *na_erk;
//	PDEAdvPlaneTS_na_sl *na_sl = nullptr;
	PDEAdvPlaneTS_interface *master;

	ShackPDEAdvectionPlaneTimeDiscretization *timeDisc;

	PDEAdvPlaneTimeSteppers()	:
		na_erk(nullptr),
		master(nullptr)
	{
	}

	void clear()
	{
		if (na_erk != nullptr)
		{
			delete na_erk;
			na_erk = nullptr;
		}
#if 0
		if (na_sl != nullptr)
		{
			delete na_sl;
			na_sl = nullptr;
		}
#endif
	}


	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		timeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}


	bool setup(
			PlaneOperators &i_op,
			sweet::ShackDictionary &io_shackDict
	)
	{
		std::cout << timeDisc->timestepping_method << std::endl;

		if (timeDisc->timestepping_method == "na_erk")
		{
			na_erk = new PDEAdvPlaneTS_na_erk(io_shackDict, i_op);
			na_erk->setup(timeDisc->timestepping_order);

			master = &(PDEAdvPlaneTS_interface&)*na_erk;
			return true;
		}
#if 0
		else if (i_timestepping_method == "na_sl")
		{
			na_sl = new Adv_Plane_TS_na_sl(i_simVars, i_op);
			na_sl->setup(i_simVars.disc.timestepping_order);

			master = &(Adv_Plane_TS_interface&)*na_sl;
			return true;
		}
#endif
		return error.set("No valid --timestepping-method provided");
	}


	~PDEAdvPlaneTimeSteppers()
	{
		clear();
	}
};




#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TIMESTEPPERS_HPP_ */
