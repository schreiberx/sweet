/*
 * SWE_Sphere_TS_ln_settls.hpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Based on plane code
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_SETTLS_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LN_SETTLS_VD_HPP_

#include <limits>
#include <string>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/core/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <sweet/core/sphere/SphereTimestepping_ExplicitRK.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"
#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"



class SWE_Sphere_TS_ln_settls_vd	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();

	std::string string_id_storage;

	void setup_auto();


private:
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &ops;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	sweet::SphereOperators_Sampler_SphereDataPhysical sphereSampler;

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

	sweet::SphereData_Spectral coriolis_arrival_spectral;
	sweet::SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;

	SWE_Sphere_TS_ln_erk_split_vd* swe_sphere_ts_ln_erk_split_vd = nullptr;
	SWE_Sphere_TS_l_irk* swe_sphere_ts_l_irk = nullptr;
	SWE_Sphere_TS_lg_irk* swe_sphere_ts_lg_irk = nullptr;



public:
	SWE_Sphere_TS_ln_settls_vd(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op,
			bool i_setup_auto = false
		);


	void setup(
			int i_timestepping_order,
			LinearCoriolisTreatment_enum i_coriolis_treatment,// = SWE_Sphere_TS_ln_settls::CORIOLIS_LINEAR,		// "ignore", "linear", "nonlinear", "semi-lagrangian"
			NLRemainderTreatment_enum i_nonlinear_divergence_treatment,// = SWE_Sphere_TS_ln_settls::NL_DIV_NONLINEAR,	// "ignore", "nonlinear"
			bool original_linear_operator_sl_treatment	// = true
	);


	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	void run_timestep_1st_order(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_ln_settls_vd();
};

#endif
