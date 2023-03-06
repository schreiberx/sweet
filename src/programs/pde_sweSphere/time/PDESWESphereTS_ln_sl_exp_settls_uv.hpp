/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_UV_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_LN_SL_EXP_SETTLS_UV_HPP_

#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/core/sphere/SphereTimestepping_SemiLagrangian.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>

#include "PDESWESphereTS_BaseInterface.hpp"
#include "PDESWESphereTS_l_exp.hpp"
#include "PDESWESphereTS_ln_erk_split_uv.hpp"



class PDESWESphereTS_ln_sl_exp_settls_uv	: public PDESWESphereTS_BaseInterface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();
	void setup_auto();

	std::string string_id_storage;

private:
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &ops;

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
	sweet::SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	sweet::SphereData_Spectral U_phi_prev, U_vrt_prev, U_div_prev;

	PDESWESphereTS_ln_erk_split_uv* swe_sphere_ts_ln_erk_split_uv = nullptr;
	PDESWESphereTS_l_exp *swe_sphere_ts_l_rexi = nullptr;


public:
	PDESWESphereTS_ln_sl_exp_settls_uv(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op,
			bool i_setup_auto = false
		);

	void print_help();

	void setup(
			int i_timestepping_order,
			LinearCoriolisTreatment_enum i_coriolis_treatment,
			NLRemainderTreatment_enum i_nonlinear_divergence_treatment,
			bool original_linear_operator_sl_treatment
	);


	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_2nd_order(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

#if (SWEET_PARAREAL && SWEET_PARAREAL_SPHERE) || (SWEET_XBRAID && SWEET_XBRAID_SPHERE)
	void set_previous_solution(
				sweet::SphereData_Spectral &i_phi_prev,
				sweet::SphereData_Spectral &i_vrt_prev,
				sweet::SphereData_Spectral &i_div_prev
	) override
	{
		if (shackDict.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		U_phi_prev = i_phi_prev;
		U_vrt_prev = i_vrt_prev;
		U_div_prev = i_div_prev;
	}
#endif

	virtual ~PDESWESphereTS_ln_sl_exp_settls_uv();
};

#endif
