/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_erk.hpp"
#include "SWE_Sphere_TS_lg_irk.hpp"



class SWE_Sphere_TS_lg_irk_lc_erk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					)
	{
		timestepping_method = i_timestepping_method;
		timestepping_order = shackDict.disc.timestepping_order;
		timestepping_order2 = shackDict.disc.timestepping_order2;
		return (
				i_timestepping_method == "lg_irk_lc_erk" ||
				i_timestepping_method == "lg_irk_lc_erk_ver0"	||
				i_timestepping_method == "lg_irk_lc_erk_ver1"
		);
	}

public:
	std::string string_id()
	{
		std::string s = "lg_irk_lc_erk_ver";

		if (version_id == 0)
			s += "0";
		else if (version_id == 1)
			s += "1";
		else
			SWEETError("Version ID");

		return s;
	}

	sweet::ShackDictionary &shackDict;
	sweet::SphereOperators &op;

	int version_id;

	int timestepping_order;
//	int timestepping_order2;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_lg_irk timestepping_lg_irk;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_lg_erk_lc_erk timestepping_lg_erk_lc_erk;

	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;


public:
	SWE_Sphere_TS_lg_irk_lc_erk(
			sweet::ShackDictionary &i_shackDict,
			sweet::SphereOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_version_id
	);


	void setup_auto();

	void run_timestep(
			sweet::SphereData_Spectral &io_phi,	///< prognostic variables
			sweet::SphereData_Spectral &io_vort,	///< prognostic variables
			sweet::SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_lg_irk_lc_erk();
};

#endif
