/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_ZERO_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_ZERO_HPP_


#include "PDEAdvectionSphere_Benchmark_BaseInterface.hpp"
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>



class PDEAdvectionSphereBenchmark_zero	:
	public PDEAdvectionSphere_Benchmark_BaseInterface
{
public:
	PDEAdvectionSphereBenchmark_zero()
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
		stream << " * ZERO FIELD TEST CASES:" << std::endl;
		stream << "    - 'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}


	void getInitialState(
		std::vector<sweet::SphereData_Spectral> &o_prognostic_fields,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		getInitialState(o_prognostic_fields[0], o_u, o_v);
	}


	void getInitialState(
		sweet::SphereData_Spectral &o_phi_pert,
		sweet::SphereData_Physical &o_u,
		sweet::SphereData_Physical &o_v
	)
	{
		o_phi_pert.spectral_set_zero();
		o_u.physical_set_zero();
		o_v.physical_set_zero();
	}
};

#endif
