/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LC_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LC_N_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_exp.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"



class SWE_Sphere_TS_lg_exp_lc_n_erk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		timestepping_order2 = shackDict.disc.timestepping_order2;
		if (
			i_timestepping_method == "lg_exp_lc_n_erk" || i_timestepping_method == "lg_exp_lc_n_erk_ver0" ||
			i_timestepping_method == "lg_exp_lc_n_erk_ver1"
		)
			return true;

		return false;
	}

	std::string string_id()
	{
		std::string s = "lg_exp_lc_n_erk_ver";

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
		if (timestepping_method == "lg_exp_lc_n_erk_ver1")
			version = 1;

		setup(
				shackDict.rexi,
				timestepping_order,
				timestepping_order2,
				shackDict.timecontrol.current_timestep_size,
				version
			);
	}


private:
	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int version_id;

	int timestepping_order;
	int timestepping_order2;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_exp timestepping_lg_rexi;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_lg_erk_lc_n_erk timestepping_lg_erk_lc_n_erk;

	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;


public:
	SWE_Sphere_TS_lg_exp_lc_n_erk(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			EXP_sweet::ShackDictionary &i_rexiSimVars,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,
			int i_version_id
	);

	void run_timestep(
			sweet::SphereData_Spectral &io_phi,		///< prognostic variables
			sweet::SphereData_Spectral &io_vrt,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_lg_exp_lc_n_erk();
};

#endif
