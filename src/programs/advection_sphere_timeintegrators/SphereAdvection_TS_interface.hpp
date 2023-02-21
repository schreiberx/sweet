/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SPHERE_ADVECTION_TS_INTERFACE_HPP_
#define SRC_PROGRAMS_SPHERE_ADVECTION_TS_INTERFACE_HPP_

#include "../advection_sphere_benchmarks/BenchmarksSphereAdvection.hpp"
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <ostream>
#include <sweet/SimulationVariables.hpp>
#include <sweet/SWEETError.hpp>


class SphereAdvection_TS_interface
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
			std::vector<SphereData_Spectral*> &io_prognostic_fields,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	)
	{
		if (io_prognostic_fields.size() == 1)
		{
			run_timestep(
					*io_prognostic_fields[0],
					io_u,
					io_v,
					i_fixed_dt,
					i_simulation_timestamp,
					i_sphereBenchmarks
				);

			return;
		}

		SWEETError("TODO: Implement time integration for multiple prognostic variables");
	}

	/*
	 * Timestepping interface used by main timestepping loop
	 */
	virtual void run_timestep(
			SphereData_Spectral &io_prognostic_field,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	)
	{
		SWEETError("TODO: Implement single prognostic variable time integration for this time integrator");
	}


	virtual bool implements_timestepping_method(
			const std::string &i_timestepping_method
		) = 0;

	virtual std::string string_id() = 0;

	virtual std::string get_help() = 0;

	virtual ~SphereAdvection_TS_interface()
	{
	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
