#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_EXP_SPECIAL_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_EXP_SPECIAL_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_lg_exp_direct.hpp"
#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"



class SWE_Sphere_TS_l_exp_direct_special	: public SWE_Sphere_TS_interface
{
	/*
	 * This class acts as a wrapper around the lg_exp_direct
	 * method and an ETDnRK version (lg + lc) to automatically
	 * choose the right one.
	 *
	 * Note, that the lc ETD method is actually not a direct method.
	 */
public:
	bool implements_timestepping_method(
		const std::string &i_timestepping_method
/////////////////#if SWEET_PARAREAL
/////////////////		,
/////////////////		int &i_timestepping_order,
/////////////////		int &i_timestepping_order2
/////////////////#endif
	)
	{
		if (i_timestepping_method == "l_exp_special")
			return true;

		if (i_timestepping_method == "lg_exp_special")
			return true;

		return false;
	}

public:

	std::string string_id()
	{
		return "l_exp_special";
	}

	SimulationVariables &simVars;
	SphereOperators_SphereData &op;
	std::string timestepping_method;

	int timestepping_order;
	bool use_coriolis;

	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_phi0;
	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_phi1;
	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_phi2;

	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_ups1;
	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_ups2;
	SWE_Sphere_TS_lg_exp_direct timestepping_lg_exp_ups3;

	SWE_Sphere_TS_ln_erk_split_vd timestepping_lc_erk;

	SWE_Sphere_TS_l_exp_direct_special(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,				///< Order of RK time stepping method
			bool i_use_coriolis,		///< Include Coriolis term
			const std::string i_function_name	///< phi/ups function
	);

	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	void euler_timestep_store_update_lc(
			const SphereData_Spectral &i_phi_pert,
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,
			SphereData_Spectral &o_phi_pert,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div,
			double i_simulation_timestamp
	);

	virtual ~SWE_Sphere_TS_l_exp_direct_special();
};

#endif

