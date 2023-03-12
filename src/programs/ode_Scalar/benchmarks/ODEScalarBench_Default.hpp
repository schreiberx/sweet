/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */
#ifndef ODE_SCALAR_BENCH_DEFAULT_BUMP_HPP__
#define ODE_SCALAR_BENCH_DEFAULT_BUMP_HPP__


#include <stdlib.h>
#include <cmath>

#include "ShackODEScalarBenchmarks.hpp"
#include "ODEScalarBench_BaseInterface.hpp"

/**
 * Setup Default benchmark
 */
class ODEScalarBenchDefault	:
		public ODEScalarBench_BaseInterface
{

public:
	bool setupBenchmark(
			double &o_u
	)
	{

		o_u = shackBenchmarks->u0;

		return true;
	}
};


#endif
