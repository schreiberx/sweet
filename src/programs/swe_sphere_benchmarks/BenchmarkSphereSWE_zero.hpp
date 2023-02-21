/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_
#define SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_


#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



class BenchmarkSphereSWE_zero	: public SWESphereBenchmarks_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	BenchmarkSphereSWE_zero()
	{

	}

	std::string benchmark_name;

	bool implements_benchmark(
			const std::string &i_benchmark_name
		)
	{
		benchmark_name = i_benchmark_name;

		return
			i_benchmark_name == "zero"
		;
	}


	void setup(
			SimulationVariables *i_simVars,
			SphereOperators_SphereData *i_ops
	)
	{
		simVars = i_simVars;
		ops = i_ops;
	}


	std::string get_help()
	{
		std::ostringstream stream;
		stream << "  'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}

	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Spectral &o_vrt,
		SphereData_Spectral &o_div
	)
	{
		o_phi_pert.spectral_set_zero();
		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
	}
};

#endif
