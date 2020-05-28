/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_NEW_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_INTERFACE_NEW_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>


class SWE_Sphere_TS_interface
{
public:

	/*
	 * Automatic setup based on simVars and operator
	 */
	virtual void setup_auto() = 0;

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void run_timestep(
			SphereData_Spectral &io_h,	///< prognostic variables
			SphereData_Spectral &io_u,	///< prognostic variables
			SphereData_Spectral &io_v,	///< prognostic variables

			double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp
	) = 0;

	virtual bool implements_timestepping_method(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string string_id() = 0;

	virtual ~SWE_Sphere_TS_interface()
	{
	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
