/*
 * Burgers_Plane_TS_ln_imex_forcing.hpp
 *
 *  Created on: 17 June 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_FORCING_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_FORCING_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <complex>
#include <sweet/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>

#include <cmath>

#include <sweet/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>

#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include <sweet/plane/PlaneStaggering.hpp>
#include <sweet/SWEETError.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"

class Burgers_Plane_TS_ln_imex_forcing	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	PlaneDataTimesteppingExplicitRK timestepping_rk;

public:
	Burgers_Plane_TS_ln_imex_forcing(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables
			PlaneData_Spectral &io_u_prev,	///< prognostic variables
			PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~Burgers_Plane_TS_ln_imex_forcing();
};

#endif /* SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_FORCING_HPP_ */
