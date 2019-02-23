/*
 * SWE_Sphere_TS_lg_erk_lf_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_N_T_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_N_T_ERK_HPP_

#include <limits>
#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereDataTimesteppingExplicitRK.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_lg_erk_lc_n_t_erk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators &op;

	int timestepping_order;
	int timestepping_order2;

	SphereDataTimesteppingExplicitRK timestepping_rk_linear;
	SphereDataTimesteppingExplicitRK timestepping_rk_nonlinear;

	// Coriolis effect
	SphereDataPhysical fg;

        // Topography
        SphereDataSpectral phi_topo;


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
	void euler_timestep_update_coriolis_and_nonlinear(
			const SphereDataSpectral &i_h,	///< prognostic variables
			const SphereDataSpectral &i_u,	///< prognostic variables
			const SphereDataSpectral &i_v,	///< prognostic variables

			SphereDataSpectral &o_h_t,	///< time updates
			SphereDataSpectral &o_u_t,	///< time updates
			SphereDataSpectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_coriolis_and_nonlinear(
			SphereDataSpectral &io_phi,		///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


private:
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
	SWE_Sphere_TS_lg_erk_lc_n_t_erk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereDataSpectral &io_phi,	///< prognostic variables
			SphereDataSpectral &io_vort,	///< prognostic variables
			SphereDataSpectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_lg_erk_lc_n_t_erk();

        const SphereDataSpectral& getTopography() const
         {
	   return phi_topo;
	 }
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
