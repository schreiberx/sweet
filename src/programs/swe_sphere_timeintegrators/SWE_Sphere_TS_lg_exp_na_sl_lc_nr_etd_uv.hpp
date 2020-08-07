/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETD_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETD_UV_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"


class SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);
	std::string string_id();
	void setup_auto();
	void print_help();

private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &ops;

	SWE_Sphere_TS_ln_erk_split_uv ts_ln_erk_split_uv;


private:
	enum NLRemainderTreatment_enum{
		NL_REMAINDER_IGNORE,
		NL_REMAINDER_NONLINEAR,
	};

	NLRemainderTreatment_enum nonlinear_remainder_treatment;

public:
	SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;

	SWE_Sphere_TS_l_exp ts_phi0_exp;
	SWE_Sphere_TS_l_exp ts_phi1_exp;
	SWE_Sphere_TS_l_exp ts_phi2_exp;

#if 0
	SWE_Sphere_TS_l_exp ts_ups0_exp;
	SWE_Sphere_TS_l_exp ts_ups1_exp;
	SWE_Sphere_TS_l_exp ts_ups2_exp;
	SWE_Sphere_TS_l_exp ts_ups3_exp;
#endif

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	int timestepping_order;
	int timestepping_order2;


public:
	SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexi,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,

			NLRemainderTreatment_enum i_nonlinear_remainder_treatment
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etd_uv();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
