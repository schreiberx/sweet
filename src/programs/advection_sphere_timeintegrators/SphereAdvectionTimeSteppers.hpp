/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_SPHERE_ADVECTION_TIMESTEPPERS_HPP_
#define SRC_SPHERE_ADVECTION_TIMESTEPPERS_HPP_

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_interface.hpp"


/**
 * SWE Plane time steppers
 */
class SphereAdvectionTimeSteppers
{
public:
	SphereAdvection_TS_interface *master = nullptr;

	std::vector<SphereAdvection_TS_interface*> registered_integrators;


	SphereAdvectionTimeSteppers();

	void reset();

	void integrators_register_all(SphereOperators_SphereData &i_op, SimulationVariables &i_simVars);

	void integrators_free_all(SphereAdvection_TS_interface *skip_this = nullptr);

	void setup(
			const std::string &i_timestepping_method,
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	);


	~SphereAdvectionTimeSteppers();
};


#endif
