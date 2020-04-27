/*
 * SWE_Sphere_TS_ln_settls.hpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_HPP_

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
#include "SWE_Sphere_TS_l_rexi.hpp"



class SWE_Sphere_TS_ln_sl_exp_settls	: public SWE_Sphere_TS_interface
{
public:
	static bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
		/*
		 * Should contain _exp and _settls
		 */
		return (
			(i_timestepping_method.find("_settls") != std::string::npos)
			&&
			(i_timestepping_method.find("_exp") != std::string::npos)
		);

		return false;
	}

	std::string string_id_storage;

	std::string string_id()
	{
		return string_id_storage;
	}


	void setup_auto()
	{
		SWE_Sphere_TS_ln_sl_exp_settls::LinearGravityTreatment_enum linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::LinearCoriolisTreatment_enum linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::NLAdvectionTreatment_enum nonlinear_advection_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_IGNORE;
		SWE_Sphere_TS_ln_sl_exp_settls::NLDivergenceTreatment_enum nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_IGNORE;

		bool original_linear_operator_sl_treatment = true;

		// Search for implicit or exp treatment of linear parts
		if (simVars.disc.timestepping_method.find("_irk_") != std::string::npos)
			linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IMPLICIT;
		else if (simVars.disc.timestepping_method.find("_exp_") != std::string::npos)
			linear_gravity_treatment = SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL;

		// Search for Coriolis
		if (simVars.disc.timestepping_method.find("l_irk") != std::string::npos || simVars.disc.timestepping_method.find("l_exp") != std::string::npos)
			linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_LINEAR;
		else if (simVars.disc.timestepping_method.find("lc_na_sl") != std::string::npos)
			linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN;
		else if (simVars.disc.timestepping_method.find("lc_") != std::string::npos)
			linear_coriolis_treatment = SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_NONLINEAR;

		if (simVars.disc.timestepping_method.find("_na_sl") != std::string::npos)
			nonlinear_advection_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN;

		// Search for Nonlinear divergence
		if (simVars.disc.timestepping_method.find("_nd_") != std::string::npos)
			nonlinear_divergence_treatment = SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_NONLINEAR;

		if (simVars.disc.timestepping_method.find("_ver2") != std::string::npos)
			original_linear_operator_sl_treatment = false;

#if 1
		string_id_storage = "";

		if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_LINEAR)
			string_id_storage += "l";
		else
			string_id_storage += "lg";

#if 0
		if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_IMPLICIT)
			string_id_storage += "_irk";
		else if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL)
			string_id_storage += "_exp";
#else
		if (linear_gravity_treatment == SWE_Sphere_TS_ln_sl_exp_settls::LINEAR_EXPONENTIAL)
			string_id_storage += "_exp";
#endif
		if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN)
			string_id_storage += "_lc";

		if (nonlinear_advection_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN)
			string_id_storage += "_na";

		if (
			linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_SEMILAGRANGIAN ||
			nonlinear_advection_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_ADV_SEMILAGRANGIAN
		)
			string_id_storage += "_sl";

		if (linear_coriolis_treatment == SWE_Sphere_TS_ln_sl_exp_settls::CORIOLIS_NONLINEAR)
			string_id_storage += "_lc";

		if (nonlinear_divergence_treatment == SWE_Sphere_TS_ln_sl_exp_settls::NL_DIV_NONLINEAR)
			string_id_storage += "_nd";

		string_id_storage += "_settls";

		if (!original_linear_operator_sl_treatment)
			string_id_storage += "_ver2";


		if (simVars.disc.timestepping_method != string_id_storage)
		{
			if (!original_linear_operator_sl_treatment)
			{
				std::cerr << "Detected time stepping method: "+string_id_storage << std::endl;
				std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
				FatalError("Autodetection of parts of time stepping methods failed!");
			}

			std::string string_id_storage2 = string_id_storage+"_ver1";
			if (simVars.disc.timestepping_method != string_id_storage2)
			{
				std::cerr << "Detected time stepping method: "+string_id_storage << std::endl;
				std::cerr << "Provided time stepping method: "+simVars.disc.timestepping_method << std::endl;
				std::cerr << "Detected alternative time stepping method: "+string_id_storage << std::endl;
				FatalError("Autodetection of parts of time stepping methods failed!");
			}
		}
#endif

		setup(
				simVars.disc.timestepping_order,
				linear_gravity_treatment,
				linear_coriolis_treatment,					// Coriolis treatment
				nonlinear_advection_treatment,		// Nonlinear advection treatment
				nonlinear_divergence_treatment,		// Nonlinear divergence treatment
				original_linear_operator_sl_treatment			// original SL linear operator treatment
			);
	}


private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

public:
	enum LinearGravityTreatment_enum {
		LINEAR_IGNORE,
		LINEAR_IMPLICIT,
		LINEAR_EXPONENTIAL,
	};


	enum LinearCoriolisTreatment_enum {
		CORIOLIS_IGNORE,
		CORIOLIS_LINEAR,
		CORIOLIS_NONLINEAR,
		CORIOLIS_SEMILAGRANGIAN,
	};


	enum NLAdvectionTreatment_enum {
		NL_ADV_IGNORE,
		NL_ADV_SEMILAGRANGIAN,
	};


	enum NLDivergenceTreatment_enum{
		NL_DIV_IGNORE,
		NL_DIV_NONLINEAR,
	};

private:
	LinearGravityTreatment_enum linear_treatment;
	LinearCoriolisTreatment_enum linear_coriolis_treatment;
	NLAdvectionTreatment_enum nonlinear_advection_treatment;
	NLDivergenceTreatment_enum nonlinear_divergence_treatment;
	int timestepping_order;
	bool original_linear_operator_sl_treatment;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical sphereSampler;

	SphereData_Spectral U_phi_prev, U_vort_prev, U_div_prev;

	SWE_Sphere_TS_l_irk* swe_sphere_ts_l_irk;
	SWE_Sphere_TS_lg_irk* swe_sphere_ts_lg_irk;
	SWE_Sphere_TS_l_rexi *swe_sphere_ts_l_rexi;


public:
	SWE_Sphere_TS_ln_sl_exp_settls(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);


	void setup(
			int i_timestepping_order,
			LinearGravityTreatment_enum i_linear_treatment,
			LinearCoriolisTreatment_enum i_coriolis_treatment,// = SWE_Sphere_TS_ln_settls::CORIOLIS_LINEAR,		// "ignore", "linear", "nonlinear", "semi-lagrangian"
			NLAdvectionTreatment_enum i_nonlinear_advection_treatment,
			NLDivergenceTreatment_enum i_nonlinear_divergence_treatment,// = SWE_Sphere_TS_ln_settls::NL_DIV_NONLINEAR,	// "ignore", "nonlinear"
			bool original_linear_operator_sl_treatment// = true
	);


	void run_timestep_pert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	void run_timestep_nonpert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


	void run_timestep_1st_order(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_ln_sl_exp_settls();
};

#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_ */
