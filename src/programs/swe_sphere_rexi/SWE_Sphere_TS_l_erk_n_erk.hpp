/*
 * SWE_Sphere_TS_ln_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_

#include <limits>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_l_erk_n_erk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators &op;

	int timestepping_order;
	int timestepping_order2;

	SphereDataTimesteppingExplicitRK timestepping_rk_linear;
	SphereDataTimesteppingExplicitRK timestepping_rk_nonlinear;

	// Coriolis effect
	SphereDataPhysical fg;


public:
	void euler_timestep_update_linear(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_nonlinear(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


private:
	void euler_timestep_update(
			const SphereData &i_phi,	///< prognostic variables
			const SphereData &i_vort,	///< prognostic variables
			const SphereData &i_div,	///< prognostic variables

			SphereData &o_phi_t,	///< time updates
			SphereData &o_vort_t,	///< time updates
			SphereData &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	SWE_Sphere_TS_l_erk_n_erk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData &io_phi,	///< prognostic variables
			SphereData &io_vort,	///< prognostic variables
			SphereData &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_erk_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
