/*
 * Adv_Sphere_TS_na_sl.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>

#include "Adv_Sphere_TS_interface.hpp"


class Adv_Sphere_TS_na_sl	: public Adv_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;

	SphereOperators_Sampler_SphereDataPhysical sampler2D;
	SphereTimestepping_SemiLagrangian semiLagrangian;

	SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;


	/**
	 * Position in physical space given in longitude/latitude angles
	 */
	ScalarDataArray posx_a, posy_a;

public:
	Adv_Sphere_TS_na_sl(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt,				///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp,

			// for varying velocity fields
			const SWESphereBenchmarksCombined *i_sphereBenchmarks
	);



	virtual ~Adv_Sphere_TS_na_sl();
};

#endif
