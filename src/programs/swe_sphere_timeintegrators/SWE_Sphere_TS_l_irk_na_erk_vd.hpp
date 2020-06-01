/*
 * SWE_Sphere_TS_l_irk_na_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_IRK_NA_ERK_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_IRK_NA_ERK_VD_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_interface.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_l_irk.hpp"
#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_erk_split_vd.hpp"



class SWE_Sphere_TS_l_irk_na_erk_vd	: public SWE_Sphere_TS_interface

{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);
	std::string string_id();

	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_irk timestepping_l_irk;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_ln_erk_split_vd timestepping_na_erk_split_vd;

	int version_id;

	int timestepping_order;
	int timestepping_order2;


public:
	SWE_Sphere_TS_l_irk_na_erk_vd(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method for linear parts
			int i_order2,	///< order of RK time stepping method for non-linear parts
			int i_version_id
	);

	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_irk_na_erk_vd();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
