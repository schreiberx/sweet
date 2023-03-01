/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_BENCHMARKS_SPHERE_ADVECTION_ZERO_HPP_
#define SRC_BENCHMARKS_SPHERE_ADVECTION_ZERO_HPP_


#include "BenchmarksSphereAdvection_interface.hpp"
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>



class BenchmarksSphereAdvection_zero	: public BenchmarksSphereAdvection_interface
{
	SimulationVariables *simVars = nullptr;
	SphereOperators_SphereData *ops = nullptr;


public:
	BenchmarksSphereAdvection_zero()
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
		stream << " * ZERO FIELD TEST CASES:" << std::endl;
		stream << "    - 'zero':	Initialize every field with 0" << std::endl;
		return stream.str();
	}


	void get_initial_state(
		std::vector<SphereData_Spectral*> &o_prognostic_fields,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		SWEETAssert(o_prognostic_fields.size() == 1, "Only scalar field supported for this benchmark!");

		get_initial_state(*o_prognostic_fields[0], o_u, o_v);
	}


	void get_initial_state(
		SphereData_Spectral &o_phi_pert,
		SphereData_Physical &o_u,
		SphereData_Physical &o_v
	)
	{
		o_phi_pert.spectral_set_zero();
		o_u.physical_set_zero();
		o_v.physical_set_zero();
	}
};

#endif
