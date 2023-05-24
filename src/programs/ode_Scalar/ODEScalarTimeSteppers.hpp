/*
 * OdeScalarTimeSteppers.hpp
 *
 *  Created on: 7th March 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_ODE_SCALAR_TIMESTEPPERS_HPP_

#include "time/ODEScalarTS_BaseInterface.hpp"
#include "time/ODEScalarTS_ln_erk.hpp"
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "time/ShackODEScalarTimeDiscretization.hpp"

/**
 * Scalar ODE time integrators
 */
class ODEScalarTimeSteppers
{
public:
	sweet::ErrorBase error;

	ODEScalarTS_ln_erk *ln_erk;
	ODEScalarTS_BaseInterface *timestepper;

	ShackODEScalarTimeDiscretization *shackTimeDisc;

	ODEScalarTimeSteppers()	:
		ln_erk(nullptr),
		timestepper(nullptr),
		shackTimeDisc(nullptr)
	{
	}

	void clear()
	{
		if (ln_erk != nullptr)
		{
			delete ln_erk;
			ln_erk = nullptr;
		}

		timestepper = nullptr;

		shackTimeDisc = nullptr;
	}


	bool shackRegistration(
			sweet::ShackDictionary &io_shackDict
	)
	{
		shackTimeDisc = io_shackDict.getAutoRegistration<ShackODEScalarTimeDiscretization>();

		return error.forwardWithPositiveReturn(io_shackDict.error);
	}


	bool setup(
			sweet::ShackDictionary &io_shackDict
	)
	{
		assert(shackTimeDisc != nullptr);

		std::cout << "Setting up time stepping method '" << shackTimeDisc->timestepping_method << "'" << std::endl;

		if (shackTimeDisc->timestepping_method == "ln_erk")
		{
			ln_erk = new ODEScalarTS_ln_erk;
			ln_erk->shackRegistration(&io_shackDict);
			ln_erk->setup();

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*ln_erk);

			timestepper = static_cast<ODEScalarTS_BaseInterface*>(ln_erk);
			return true;
		}

		return error.set("No valid --timestepping-method=... provided ('"+shackTimeDisc->timestepping_method+"')");
	}


	~ODEScalarTimeSteppers()
	{
		clear();
	}
};




#endif
