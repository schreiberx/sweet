/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../advection_sphere_timeintegrators/SphereAdvectionTimeSteppers.hpp"

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_erk.hpp"
#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_sl.hpp"
#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_trajectories.hpp"



SphereAdvectionTimeSteppers::SphereAdvectionTimeSteppers()
{
}

void SphereAdvectionTimeSteppers::reset()
{
	delete master;
	master = nullptr;

	integrators_free_all();
}



void SphereAdvectionTimeSteppers::integrators_register_all(SphereOperators_SphereData &i_op, SimulationVariables &i_simVars)
{
	/*
	 * Register time integrators
	 */
	registered_integrators.push_back(static_cast<SphereAdvection_TS_interface*>(new SphereAdvection_TS_na_erk(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SphereAdvection_TS_interface*>(new SphereAdvection_TS_na_sl(i_simVars, i_op)));
	registered_integrators.push_back(static_cast<SphereAdvection_TS_interface*>(new SphereAdvection_TS_na_trajectories(i_simVars, i_op)));
}



void SphereAdvectionTimeSteppers::integrators_free_all(SphereAdvection_TS_interface *skip_this)
{

	for (std::size_t i = 0; i < registered_integrators.size(); i++)
	{
		SphereAdvection_TS_interface *ts = registered_integrators[i];

		if (ts == skip_this)
			continue;

		delete ts;
	}

	registered_integrators.clear();
}



void SphereAdvectionTimeSteppers::setup(const std::string &i_timestepping_method, SphereOperators_SphereData &i_op, SimulationVariables &i_simVars)
{
	reset();

	integrators_register_all(i_op, i_simVars);

	/*
	 * Find right one
	 */
	master = nullptr;

	for (std::size_t i = 0; i < registered_integrators.size(); i++)
	{
		SphereAdvection_TS_interface *ts = registered_integrators[i];

		if (ts->implements_timestepping_method(i_timestepping_method))
		{
			if (master != nullptr)
			{
				std::cout << "Processing " << i+1 << "th element" << std::endl;
				SWEETError(std::string("Duplicate implementation for method ") + i_timestepping_method);
			}

			std::cout << "Found match at " << i+1 << "th element" << std::endl;
			ts->setup_auto();
			master = ts;
		}
	}

	if (master == nullptr)
		SWEETError(std::string("No valid --timestepping-method '") + i_timestepping_method + std::string("' provided"));

	// Found integrator, freeing others
	integrators_free_all(master);
}


SphereAdvectionTimeSteppers::~SphereAdvectionTimeSteppers()
{
	reset();
}
