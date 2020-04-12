/*
 * SWE_Sphere_TS_l_rexi_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SWE_SPHERE_TS_L_REXI_N_ERK_HPP_
#define SWE_SPHERE_TS_L_REXI_N_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"



class SWE_Sphere_TS_l_rexi_n_erk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	double timestep_size;

	/*
	 * Linear time steppers
	 */
	SWE_Sphere_TS_l_rexi timestepping_l_rexi;

	/*
	 * Non-linear time steppers
	 */
	SWE_Sphere_TS_l_erk_n_erk timestepping_l_erk_n_erk;

	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;

	int version_id;

	int timestepping_order;
	int timestepping_order2;


public:
	SWE_Sphere_TS_l_rexi_n_erk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexiSimVars,
			int i_order,	///< order of RK time stepping method
			int i_order2,	///< order of RK time stepping method
			double i_timestep_size,
			bool i_use_f_sphere,
			int version_id	///< strang splitting for 2nd order method
	);

	void run_timestep_pert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	void run_timestep_nonpert(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_rexi_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
