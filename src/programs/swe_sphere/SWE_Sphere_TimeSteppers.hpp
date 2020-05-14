/*
 * SWE_Sphere_TimeSteppers.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
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


	SWE_Sphere_TimeSteppers();

	void reset();

	void setup(
			const std::string &i_timestepping_method,
			SphereOperators_SphereData &i_op,
			SimulationVariables &i_simVars
	);


	~SWE_Sphere_TimeSteppers();
};




#endif
