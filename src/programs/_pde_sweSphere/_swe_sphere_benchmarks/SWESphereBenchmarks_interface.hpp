/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_
#define SRC_SWE_SPHERE_BENCHMARKS_TS_INTERFACE_HPP_

#include <sweet/core/ScalarDataArray.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>


class SWESphereBenchmarks_interface
{
public:
	virtual void setup(
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *i_ops
		) = 0;

	virtual bool implements_benchmark(
			const std::string &i_benchmark_name
		) = 0;

	virtual void get_initial_state(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
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
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div,
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

#endif
