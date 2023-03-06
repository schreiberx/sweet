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

void PDEAdvectionSphereTimeSteppers::setup_1_registerAllTimesteppers()
{
	/*
	 * Register time integrators
	 */
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_erk));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_sl));
	_registered_integrators.push_back(static_cast<PDEAdvectionSphereTS_BaseInterface*>(new PDEAdvectionSphereTS_na_trajectories));
}



bool PDEAdvectionSphereTimeSteppers::setup_2_shackRegistration(
		sweet::ShackDictionary *io_shackDict
)
{
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->shackRegistration(io_shackDict);
	}
	return true;
}



void PDEAdvectionSphereTimeSteppers::printImplementedTimesteppingMethods(
	std::ostream &o_ostream,
	const std::string &i_prefix
)
{
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (START)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;

	std::string prefix = i_prefix+"  ";
	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		_registered_integrators[i]->printImplementedTimesteppingMethods(o_ostream, prefix);
	}

	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << "Timestepping methods (END)" << std::endl;
	o_ostream << "********************************************************************************" << std::endl;
	o_ostream << std::endl;
}


bool PDEAdvectionSphereTimeSteppers::setup_3_timestepper(
		const std::string &i_timestepping_method,
		sweet::ShackDictionary *io_shackDict,
		sweet::SphereOperators *io_ops
)
{
	if (i_timestepping_method == "")
	{
		printImplementedTimesteppingMethods();
		return error.set("Please set time stepping method using --timestepping-method=...");
	}
	/*
	 * Find right one
	 */
	timestepper = nullptr;

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts->testImplementsTimesteppingMethod(i_timestepping_method))
		{
			if (timestepper != nullptr)
			{
				//std::cout << "Processing " << i+1 << "th element" << std::endl;
				return error.set("Duplicate implementation for method "+i_timestepping_method);
			}

			//std::cout << "Found matching time stepping method at " << i+1 << "th element" << std::endl;
			ts->setup(io_ops);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(*ts);
			timestepper = ts;
		}
	}

	if (timestepper == nullptr)
		return error.set("No valid --timestepping-method '"+i_timestepping_method+"' provided");

	// Found integrator, freeing others
	_timesteppersFreeAll(timestepper);

	return true;
}


void PDEAdvectionSphereTimeSteppers::_timesteppersFreeAll(
		PDEAdvectionSphereTS_BaseInterface *i_skip_this_timestepper
)
{

	for (std::size_t i = 0; i < _registered_integrators.size(); i++)
	{
		PDEAdvectionSphereTS_BaseInterface *ts = _registered_integrators[i];

		if (ts == i_skip_this_timestepper)
			continue;

		delete ts;
	}

	_registered_integrators.clear();
}


void PDEAdvectionSphereTimeSteppers::clear()
{
	delete timestepper;
	timestepper = nullptr;

	_timesteppersFreeAll();
}


PDEAdvectionSphereTimeSteppers::~PDEAdvectionSphereTimeSteppers()
{
	clear();
}
