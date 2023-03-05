/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#include "time/PDEAdvectionSphereTS_na_erk.hpp"
#include "time/PDEAdvectionSphereTS_na_sl.hpp"
#include "time/PDEAdvectionSphereTS_na_trajectories.hpp"
#include "PDEAdvectionSphereTimeSteppers.hpp"




PDEAdvectionSphereTimeSteppers::PDEAdvectionSphereTimeSteppers()
{
}


void PDEAdvectionSphereTimeSteppers::_integratorsRegisterAll(
		sweet::ShackDictionary *i_shackDict,
		sweet::SphereOperators *i_op
)
{
	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_erk));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_sl));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_trajectories));
}



void PDEAdvectionSphereTimeSteppers::_integratorsFreeAll(PDEAdvectionSphereTS_BaseInterface *skip_this)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}



bool PDEAdvectionSphereTimeSteppers::setup(
		const std::string &i_timestepping_method,
		sweet::ShackDictionary *io_shackDict,
		sweet::SphereOperators *io_ops
)
{
	_integratorsRegisterAll(io_shackDict, io_ops);

	/*
	 * Find right one
	 */
	master = nullptr;

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->implements_timestepping_method(i_timestepping_method))
		{
			if (master != nullptr)
			{
				std::cout << "Processing " << i+1 << "th element" << std::endl;
				SWEETError(std::string("Duplicate implementation for method ") + i_timestepping_method);
			}

			std::cout << "Found matching time stepping method at " << i+1 << "th element" << std::endl;
			ts->setup(io_ops);
			master = ts;
		}
	}

	if (master == nullptr)
		SWEETError(std::string("No valid --timestepping-method '") + i_timestepping_method + std::string("' provided"));

	// Found integrator, freeing others
	_integratorsFreeAll(master);

	return true;
}

void PDEAdvectionSphereTimeSteppers::clear()
{
	delete master;
	master = nullptr;

	_integratorsFreeAll();
}

bool PDEAdvectionSphereTimeSteppers::shackRegistration(
	sweet::ShackDictionary *io_shackDict
)
{
	PDEAdvectionSphereTS_na_erk a;
	a.shackRegistration(io_shackDict);

	ERROR_CHECK_WITH_RETURN_BOOLEAN(a);
	return true;
}


PDEAdvectionSphereTimeSteppers::~PDEAdvectionSphereTimeSteppers()
{
	clear();
}
