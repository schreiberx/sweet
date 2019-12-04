/*
 * SWE_bench_NormalModes.hpp
 *
 *  Created on: 03 Nov 2019
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 */
#ifndef SWE_PLANE_NORMAL_MODES_HPP_
#define SWE_PLANE_NORMAL_MODES_HPP_

#include <stdlib.h>

#include <sweet/sweetmath.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <libmath/GaussQuadrature.hpp>


/**
 * Implement normal mode initialization
 *
 *
 *
 **/
class SWE_bench_NormalModes
{
	SimulationVariables &simVars;

	PlaneOperators &op;

	double f = simVars.sim.plane_rotating_f0;
	double g = simVars.sim.gravitation;
	double sx = simVars.sim.plane_domain_size[0];
	double sy = simVars.sim.plane_domain_size[1];



public:
	SWE_bench_NormalModes(
		SimulationVariables &io_simVars,
		PlaneOperators &io_op
	)	:
		simVars(io_simVars),
		op(io_op)
	{
	}

	void setup(
			PlaneData &o_h,
			PlaneData &o_u,
			PlaneData &o_v
	)
	{
		std::cout<< "Generating Normal Modes Initial Conditions: ";

		if(simVars.benchmark.benchmark_normal_modes_case==""){
			FatalError("SWE_bench_NormalModes: please choose the normal mode case with --benchmark-normal-mode-case [string] (see SWE_bench_NormalModes.hpp file)");
		};
		std::cout<< simVars.benchmark.benchmark_normal_modes_case;

		if(simVars.benchmark.benchmark_normal_modes_case=="single")
		{

		};
		
		std::cout<< "   Done! " << std::endl;
	}


};


#endif
