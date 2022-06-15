/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "SWE_Sphere_TS_interface.hpp"

/**
 * SWE Plane time steppers
 */
class SWE_Sphere_TimeSteppers
{
public:
	SWE_Sphere_TS_interface *master = nullptr;

	std::vector<SWE_Sphere_TS_interface*> registered_integrators;


	SWE_Sphere_TimeSteppers();

	void reset();

	void integrators_register_all(SphereOperators_SphereData &i_op, SimulationVariables &i_simVars);

	void integrators_free_all(SWE_Sphere_TS_interface *skip_this = nullptr);

	void print_help_all_registered();

	void setup(
			const std::string &i_timestepping_method,
#if SWEET_PARAREAL
			int &i_timestepping_order,
			int &i_timestepping_order2,
#endif
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	);

	~SWE_Sphere_TimeSteppers();
};




#endif
