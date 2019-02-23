/*
 * SWE_Sphere_TS_ln_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_L_ERK_N_ERK_HPP_

#include <limits>
#include <sweet/sphere/SphereDataSpectral.hpp>
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
			const SphereDataSpectral &i_h,	///< prognostic variables
			const SphereDataSpectral &i_u,	///< prognostic variables
			const SphereDataSpectral &i_v,	///< prognostic variables

			SphereDataSpectral &o_h_t,	///< time updates
			SphereDataSpectral &o_u_t,	///< time updates
			SphereDataSpectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_n(
			const SphereDataSpectral &i_h,	///< prognostic variables
			const SphereDataSpectral &i_u,	///< prognostic variables
			const SphereDataSpectral &i_v,	///< prognostic variables

			SphereDataSpectral &o_h_t,	///< time updates
			SphereDataSpectral &o_u_t,	///< time updates
			SphereDataSpectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_n(
			SphereDataSpectral &io_phi,	///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

public:
	void euler_timestep_update(
			const SphereDataSpectral &i_phi,	///< prognostic variables
			const SphereDataSpectral &i_vort,	///< prognostic variables
			const SphereDataSpectral &i_div,	///< prognostic variables

			SphereDataSpectral &o_phi_t,	///< time updates
			SphereDataSpectral &o_vort_t,	///< time updates
			SphereDataSpectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	SWE_Sphere_TS_l_erk_n_erk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2
	);

	void run_timestep(
			SphereDataSpectral &io_phi,	///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_l_erk_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
