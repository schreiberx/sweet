/*
 * Burgers_Plane_TS_ln_imex.hpp
 *
 *  Created on: 17 June 2017
 *  Author: Andreas Schmitt <aschmitt@fnb.tu-darmstadt.de>
 */

#ifndef SRC_PROGRAMS_Burgers_PLANE_TS_LN_IMEX_HPP_
#define SRC_PROGRAMS_Burgers_PLANE_TS_LN_IMEX_HPP_

#include <limits>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include "Burgers_Plane_TS_interface.hpp"

#include <complex>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperatorsComplex.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>

#include <sweet/sweetmath.hpp>

#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>

#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>
#include <sweet/FatalError.hpp>

#include <sweet/plane/Staggering.hpp>

class Burgers_Plane_TS_ln_imex	: public Burgers_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	PlaneDataTimesteppingRK timestepping_rk;

public:
	Burgers_Plane_TS_ln_imex(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables
			PlaneData &io_u_prev,	///< prognostic variables
			PlaneData &io_v_prev,	///< prognostic variables

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1,
			double i_max_simulation_time = std::numeric_limits<double>::infinity()
	);



	virtual ~Burgers_Plane_TS_ln_imex();
};

#endif /* SRC_PROGRAMS_Burgers_PLANE_TS_LN_IMEX_HPP_ */
