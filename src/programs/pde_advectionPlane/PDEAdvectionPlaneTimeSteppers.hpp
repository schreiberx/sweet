/*
 * Adv_Plane_TimeSteppers.hpp
 *
 *  Created on: 4th April 2018
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTION_PLANE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_PDE_ADVECTION_PLANE_TIMESTEPPERS_HPP_

#include "time/PDEAdvectionPlaneTS_BaseInterface.hpp"
#include "time/PDEAdvectionPlaneTS_na_erk.hpp"
#include "time/PDEAdvectionPlaneTS_na_sl.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "time/ShackPDEAdvectionPlaneTimeDiscretization.hpp"

/**
 * Advection for the plane time integrators
 */
class PDEAdvectionPlaneTimeSteppers
{
public:
	sweet::ErrorBase error;

	PDEAdvectionPlaneTS_na_erk *na_erk;
	PDEAdvectionPlaneTS_na_sl *na_sl;
	PDEAdvectionPlaneTS_BaseInterface *master;

	ShackPDEAdvectionPlaneTimeDiscretization *shackTimeDisc;

	PDEAdvectionPlaneTimeSteppers()	:
		na_erk(nullptr),
		na_sl(nullptr),
		master(nullptr),
		shackTimeDisc(nullptr)
	{
	}

	void clear()
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

		master = nullptr;

		shackTimeDisc = nullptr;
	}


	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackTimeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}


	bool setup(
			sweet::ShackDictionary &io_shackDict,
			sweet::PlaneOperators &i_op
	)
	{
		assert(shackTimeDisc != nullptr);

		std::cout << "Setting up time stepping method '" << shackTimeDisc->timestepping_method << "'" << std::endl;

		if (shackTimeDisc->timestepping_method == "na_erk")
		{
			na_erk = new PDEAdvectionPlaneTS_na_erk;
			na_erk->shackRegistration(&io_shackDict);
			na_erk->setup(&i_op);

			ERROR_CHECK_WITH_RETURN_BOOLEAN(*na_erk);

			master = static_cast<PDEAdvectionPlaneTS_BaseInterface*>(na_erk);
			return true;
		}
		else if (shackTimeDisc->timestepping_method == "na_sl")
		{
			na_sl = new PDEAdvectionPlaneTS_na_sl;
			na_sl->shackRegistration(&io_shackDict);
			na_sl->setup(&i_op);

			ERROR_CHECK_WITH_RETURN_BOOLEAN(*na_sl);

			master = static_cast<PDEAdvectionPlaneTS_BaseInterface*>(na_sl);
			return true;
		}

		return error.set("No valid --timestepping-method=... provided ('"+shackTimeDisc->timestepping_method+"')");
	}


	~PDEAdvectionPlaneTimeSteppers()
	{
		clear();
	}
};




#endif
