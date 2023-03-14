/*
 * Burgers_Plane_TS_ln_imex_forcing.hpp
 *
 *  Created on: 17 June 2017
 * Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_FORCING_HPP_
#define SRC_PROGRAMS_BURGERS_PLANE_TS_LN_IMEX_FORCING_HPP_

#include <limits>
#include <sweet/core/plane/sweet::PlaneData_Spectral.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <complex>
#include <sweet/core/plane/sweet::PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>

#include <cmath>

#include <sweet/core/plane/Convert_PlaneDataSpectral_to_PlaneDataSpectralComplex.hpp>
#include <sweet/core/plane/Convert_PlaneDataSpectralComplex_to_PlaneDataSpectral.hpp>

#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include <sweet/core/plane/PlaneStaggering.hpp>
#include <sweet/core/SWEETError.hpp>
#include "Burgers_Plane_TS_interface.hpp"
#include "../burgers_benchmarks/BurgersValidationBenchmarks.hpp"

class Burgers_Plane_TS_ln_imex_forcing	: public Burgers_Plane_TS_interface
{
	sweet::ShackDictionary &shackDict;
	PlaneOperators &op;

	int timestepping_order;
	PlaneDataTimesteppingExplicitRK timestepping_rk;

public:
	Burgers_Plane_TS_ln_imex_forcing(
			sweet::ShackDictionary &i_shackDict,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void runTimestep(
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_u_prev,	///< prognostic variables
			///sweet::PlaneData_Spectral &io_v_prev,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~Burgers_Plane_TS_ln_imex_forcing();
};

#endif
