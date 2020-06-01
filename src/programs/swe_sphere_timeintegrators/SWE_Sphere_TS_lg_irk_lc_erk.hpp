/*
 * SWE_Sphere_TS_lg_irk_lf_erk.hpp
 *
 *  Created on: 11 Nov 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_IRK_LC_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_interface.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_erk_lc_erk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_lg_irk.hpp"



class SWE_Sphere_TS_lg_irk_lc_erk	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method)
	{
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

	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

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
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_version_id
	);


	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Sphere_TS_lg_irk_lc_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
