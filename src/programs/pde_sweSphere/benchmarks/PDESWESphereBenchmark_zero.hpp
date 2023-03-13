/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_
#define SRC_BENCHMARK_SPHERE_SWE_ZERO_HPP_


#include "PDESWESphereBenchmarks_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class PDESWESphereBenchmark_zero	: public PDESWESphereBenchmarks_BaseInterface
{
public:
	PDESWESphereBenchmark_zero()
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


	void setup_1_shackData()
	{
	}


	void setup_2_withOps(
			sweet::SphereOperators *io_ops
	)
	{
		ops = io_ops;
	}


	void clear()
	{
	}


	std::string printHelp()
	{
		std::ostringstream stream;
		stream << "  'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}

	void getInitialState(
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
