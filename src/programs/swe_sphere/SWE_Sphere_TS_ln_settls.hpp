/*
 * SWE_Sphere_TS_ln_settls.hpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Babed on plane code
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_lg_erk.hpp"



class SWE_Sphere_TS_ln_settls	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	bool include_nonlinear_divergence;
	bool original_linear_operator_sl_tretment;
	int coriolis_treatment;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical sphereSampler;

	SphereData_Spectral phi_prev, vort_prev, div_prev;

	// Arrival points for semi-lag
	ScalarDataArray pos_lon_a, pos_lat_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	SWE_Sphere_TS_l_erk* swe_sphere_ts_l_erk;
	SWE_Sphere_TS_lg_erk* swe_sphere_ts_lg_erk;
	SWE_Sphere_TS_l_irk* swe_sphere_ts_l_irk;
	SWE_Sphere_TS_lg_irk* swe_sphere_ts_lg_irk;

	// Coriolis effect
	SphereData_Physical fg;


public:
	SWE_Sphere_TS_ln_settls(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			bool include_nonlinear_divergence = true,
			bool original_linear_operator_sl_tretment = true,
			const std::string &i_coriolis_treatment = "nonlinear"	// "linear", "nonlinear", "semi-lagrangian"
	);

	enum{
		CORIOLIS_LINEAR,
		CORIOLIS_NONLINEAR,
		CORIOLIS_SEMILAGRANGIAN,
	};

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_ln_settls();
};

#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_ */
