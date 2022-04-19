/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPERS_HPP_

#include "SWE_Sphere_TS_interface.hpp"

#if SWEET_PARAREAL
#include <sweet/plane/PlaneOperators.hpp>
#endif

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
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	);

#if SWEET_PARAREAL
	void setup(
			const std::string &i_timestepping_method,
			int &i_timestepping_order,
			int &i_timestepping_order2,

			PlaneOperators &i_op_plane,
			SphereOperators_SphereData &i_op_sphere,
			SimulationVariables &i_simVars
	);
#endif


	~SWE_Sphere_TimeSteppers();
};




#endif
