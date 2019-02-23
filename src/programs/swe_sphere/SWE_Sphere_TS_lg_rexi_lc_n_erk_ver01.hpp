/*
 * SWE_Sphere_TS_lg_rexi_lf_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LF_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_REXI_LF_N_ERK_HPP_

#include <limits>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"



class SWE_Sphere_TS_lg_rexi_lc_n_erk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators &op;

	int version_id;

	int timestepping_order;
	int timestepping_order2;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_rexi timestepping_lg_rexi;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_lg_erk_lc_n_erk timestepping_lg_erk_lc_n_erk;

	SphereDataTimesteppingExplicitRK timestepping_rk_nonlinear;


public:
	SWE_Sphere_TS_lg_rexi_lc_n_erk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexiSimVars,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size,
			int i_version_id
	);

	void run_timestep(
			SphereDataSpectral &io_phi,		///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,		///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_lg_rexi_lc_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
