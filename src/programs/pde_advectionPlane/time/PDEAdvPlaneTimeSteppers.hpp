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
#include "PDEAdvPlaneTS_na_sl.hpp"

#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "ShackPDEAdvectionPlaneTimeDiscretization.hpp"

/**
 * Advection for the plane time integrators
 */
class PDEAdvPlaneTimeSteppers
{
public:
	sweet::ErrorBase error;

	PDEAdvPlaneTS_na_erk *na_erk;
	PDEAdvPlaneTS_na_sl *na_sl;
	PDEAdvPlaneTS_interface *master;

	ShackPDEAdvectionPlaneTimeDiscretization *shackTimeDisc;

	PDEAdvPlaneTimeSteppers()	:
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
		// Get this shack back again in case a clear() followed a setup() without the shackRegistration()
		shackTimeDisc = io_shackDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();

		std::cout << "Setting up time stepping method '" << shackTimeDisc << "'" << std::endl;
		std::cout << "Setting up time stepping method '" << shackTimeDisc->timestepping_method << "'" << std::endl;

		if (shackTimeDisc->timestepping_method == "na_erk")
		{
			na_erk = new PDEAdvPlaneTS_na_erk(io_shackDict, i_op);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(*na_erk);

			master = &(PDEAdvPlaneTS_interface&)*na_erk;
			return true;
		}
		else if (shackTimeDisc->timestepping_method == "na_sl")
		{
			na_sl = new PDEAdvPlaneTS_na_sl(io_shackDict, i_op);

			master = &(PDEAdvPlaneTS_interface&)*na_sl;
			return true;
		}

		return error.set("No valid --timestepping-method=... provided ('"+shackTimeDisc->timestepping_method+"')");
	}


	~PDEAdvPlaneTimeSteppers()
	{
		clear();
	}
};




#endif
