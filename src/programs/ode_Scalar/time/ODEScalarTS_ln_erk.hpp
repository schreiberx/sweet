/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PDE_ODE_SCALAR_TS_LN_ERK_HPP_
#define SRC_PDE_ODE_SCALAR_TS_LN_ERK_HPP_

#include "ODEScalarTS_BaseInterface.hpp"
#include <sweet/core/ErrorBase.hpp>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "ShackODEScalarTimeDiscretization.hpp"
#include "../benchmarks/ShackODEScalarBenchmarks.hpp"

class ODEScalarTS_ln_erk	:
		public ODEScalarTS_BaseInterface
{
public:
	double param_a;
	double param_b;

public:
	ODEScalarTS_ln_erk();

	~ODEScalarTS_ln_erk();

	bool setup();


public:
	void runTimestep(
			double &io_u,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);
};

#endif
