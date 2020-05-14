/*
 * SWE_Sphere_TS_ln_settls.hpp
 *
 *  Created on: 24 Sep 2019
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_UV_HPP_

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
#include "SWE_Sphere_TS_ln_erk_split_uv.hpp"



class SWE_Sphere_TS_ln_sl_exp_settls_uv	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
		/*
		 * Should contain _exp and _settls as well as _uv to indicate vorticity-divergence formulation
		 */
		return (
			(i_timestepping_method.find("_settls") != std::string::npos)
			&&
			(i_timestepping_method.find("_uv") != std::string::npos)
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


	void setup_auto();

private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

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

	SWE_Sphere_TS_ln_erk_split_uv* swe_sphere_ts_ln_erk_split_uv = nullptr;
	SWE_Sphere_TS_l_exp *swe_sphere_ts_l_rexi = nullptr;


public:
	SWE_Sphere_TS_ln_sl_exp_settls_uv(
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


	void interpolate_departure_point_uv(
			const SphereData_Spectral &i_phi,
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,

			const ScalarDataArray &i_pos_lon_d,
			const ScalarDataArray &i_pos_lat_d,

			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	);

	void run_timestep_pert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_ln_sl_exp_settls_uv();
};

#endif /* SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_CN_NA_SL_ND_SETTLS_HPP_ */
