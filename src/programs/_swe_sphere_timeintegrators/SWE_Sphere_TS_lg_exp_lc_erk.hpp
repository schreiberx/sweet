/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPER_LG_EXP_LC_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TIMESTEPPER_LG_EXP_LC_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <sweet/core/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/core/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"



class SWE_Sphere_TS_lg_exp_lc_erk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = simVars.disc.timestepping_order;
		timestepping_order2 = simVars.disc.timestepping_order2;
		if (
			i_timestepping_method == "lg_exp_lc_erk" || i_timestepping_method == "lg_exp_lc_erk_ver0" ||
			i_timestepping_method == "lg_exp_lc_erk_ver1"
		)
			return true;

		return false;
	}

	std::string string_id()
	{
		std::string s = "lg_exp_lc_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETError("Version ID");

		return s;
	}

	void setup_auto()
	{
		int version = 0;
		if (timestepping_method == "lg_exp_lc_erk_ver1")
			version = 1;

		setup(
				simVars.rexi,
				timestepping_order,
				timestepping_order2,
				simVars.timecontrol.current_timestep_size,
				version
			);
	}


private:
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int version_id;

	int timestepping_order;
	int timestepping_order2;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_exp timestepping_l_exp;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_lg_erk_lc_erk timestepping_lg_erk_lc_erk;

	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;


public:
	SWE_Sphere_TS_lg_exp_lc_erk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			EXP_SimulationVariables &i_rexiSimVars,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,
			int i_version_id
	);

	void run_timestep(
			SphereData_Spectral &io_phi,		///< prognostic variables
			SphereData_Spectral &io_vrt,	///< prognostic variables
			SphereData_Spectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_lg_exp_lc_erk();
};

#endif
