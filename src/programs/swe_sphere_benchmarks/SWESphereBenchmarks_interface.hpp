/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_

#include <sweet/ScalarDataArray.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>


class SWESphereBenchmarks_interface
{
public:
	virtual void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
		) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;

	virtual void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	) = 0;

	virtual std::string get_help() = 0;

	/*
	 * Setup topography in Simulation Variables class
	 */
	virtual void setup_topography()
	{
		SWEETError("Not implemented for this benchmark");
	}

	virtual void get_reference_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div,
		double i_timestamp
	)
	{
		SWEETError("Not implemented for this benchmark");
	}

	virtual bool has_time_varying_state()
	{
		return false;
	}

	virtual ~SWESphereBenchmarks_interface()
	{
	}
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
