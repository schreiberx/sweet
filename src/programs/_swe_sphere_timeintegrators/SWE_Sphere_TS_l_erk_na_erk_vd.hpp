/*
 * SWE_Sphere_TS_l_erk_na_erk_vd.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_VD_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_L_ERK_NA_ERK_VD_HPP_

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/core/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_ln_erk_split_vd.hpp"



class SWE_Sphere_TS_l_erk_na_erk_vd	: public SWE_Sphere_TS_interface
{
public:
	bool implements_timestepping_method(const std::string &i_timestepping_method
					);
	std::string string_id();


	SimulationVariables &simVars;
	SphereOperators_SphereData &ops;

	int timestepping_order;
	int timestepping_order2;

	SWE_Sphere_TS_ln_erk_split_vd *l_erk_split_vd = nullptr;
	SWE_Sphere_TS_ln_erk_split_vd *na_erk_split_vd = nullptr;

public:
	SWE_Sphere_TS_l_erk_na_erk_vd(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2
	);


	void setup_auto();

	void run_timestep(
			SphereData_Spectral &io_phi,
			SphereData_Spectral &io_vrt,
			SphereData_Spectral &io_div,

			double i_fixed_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Sphere_TS_l_erk_na_erk_vd();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
