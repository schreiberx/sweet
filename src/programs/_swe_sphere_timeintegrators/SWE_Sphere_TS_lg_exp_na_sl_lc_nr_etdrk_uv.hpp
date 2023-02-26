/*
 * SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv
 *
 * Created on: 24 Mar 2022
 *     Author: Joao Steinstraesser <joao.steinstraesser@usp.br>
 *
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETDRK_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LG_EXP_NA_SL_LC_N_ETDRK_UV_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/core/SimulationVariables.hpp>

#include <sweet/core/sphere/SphereTimestepping_SemiLagrangian.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"


class SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv	: public SWE_Sphere_TS_interface
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
	SWE_Sphere_TS_l_exp ts_phi2_exp;

	SWE_Sphere_TS_l_exp ts_psi1_exp;
	SWE_Sphere_TS_l_exp ts_psi2_exp;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	int timestepping_order;
	int timestepping_order2;


public:
	SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			EXP_SimulationVariables &i_rexi,
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

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE)
	void set_previous_solution(
				SphereData_Spectral &i_phi_prev,
				SphereData_Spectral &i_vrt_prev,
				SphereData_Spectral &i_div_prev
	) override
	{
		if (simVars.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		U_phi_prev = i_phi_prev;
		U_vrt_prev = i_vrt_prev;
		U_div_prev = i_div_prev;
	}
#endif

	virtual ~SWE_Sphere_TS_lg_exp_na_sl_lc_nr_etdrk_uv();
};

#endif
