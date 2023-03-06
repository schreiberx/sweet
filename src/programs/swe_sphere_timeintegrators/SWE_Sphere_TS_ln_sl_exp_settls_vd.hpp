/*
 * SWE_Sphere_TS_ln_settls.hpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_VD_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"



class SWE_Sphere_TS_ln_sl_exp_settls_vd	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();
	void setup_auto();
	void print_help();

	std::string string_id_storage;

private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &ops;

public:
	enum LinearCoriolisTreatment_enum {
		CORIOLIS_IGNORE,
		CORIOLIS_LINEAR,
		CORIOLIS_NONLINEAR,
		CORIOLIS_SEMILAGRANGIAN,
	};

	enum NLRemainderTreatment_enum{
		NL_REMAINDER_IGNORE,
		NL_REMAINDER_NONLINEAR,
	};

private:
	LinearCoriolisTreatment_enum coriolis_treatment;
	NLRemainderTreatment_enum nonlinear_remainder_treatment;

	int timestepping_order;
	bool original_linear_operator_sl_treatment;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;

	SWE_Sphere_TS_ln_erk_split_vd* swe_sphere_ts_ln_erk_split_vd = nullptr;
	SWE_Sphere_TS_l_exp *swe_sphere_ts_l_exp = nullptr;


public:
	SWE_Sphere_TS_ln_sl_exp_settls_vd(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op,
			bool i_setup_auto = false
		);


	void setup(
			int i_timestepping_order,
			LinearCoriolisTreatment_enum i_coriolis_treatment,
			NLRemainderTreatment_enum i_nonlinear_divergence_treatment,
			bool original_linear_operator_sl_treatment
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
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

	virtual ~SWE_Sphere_TS_ln_sl_exp_settls_vd();
};

#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_ */
