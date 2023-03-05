/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_
#define SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_


#include "SWESphereBenchmarks_interface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class BenchmarkSphereSWE_zero	: public SWESphereBenchmarks_interface
{
	sweet::ShackDictionary *shackDict = nullptr;
	sweet::SphereOperators *ops = nullptr;


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
			sweet::ShackDictionary *i_shackDict,
			sweet::SphereOperators *i_ops
	)
	{
		shackDict = i_shackDict;
		ops = i_ops;
	}


	std::string get_help()
	{
		std::ostringstream stream;
		stream << "  'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}

	void get_initial_state(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Spectral &o_vrt,
		sweet::SphereData_Spectral &o_div
	)
	{
		o_phi_pert.spectral_set_zero();
		o_vrt.spectral_set_zero();
		o_div.spectral_set_zero();
	}
};

#endif
