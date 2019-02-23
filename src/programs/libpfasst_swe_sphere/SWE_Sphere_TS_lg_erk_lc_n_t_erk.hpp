/*
 * SWE_Sphere_TS_lg_erk_lf_n_erk.hpp
 *
 *  Created on: 30 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_N_T_ERK_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_REXI_SWE_SPHERE_TS_LG_ERK_LC_N_T_ERK_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_ExplicitRK.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include "SWE_Sphere_TS_interface.hpp"



class SWE_Sphere_TS_lg_erk_lc_n_t_erk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;
	int timestepping_order2;

	SphereTimestepping_ExplicitRK timestepping_rk_linear;
	SphereTimestepping_ExplicitRK timestepping_rk_nonlinear;

	// Coriolis effect
	SphereData_Physical fg;

        // Topography
        SphereData_Spectral phi_topo;


public:
	void euler_timestep_update_linear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_coriolis_and_nonlinear(
			const SphereData_Spectral &i_h,	///< prognostic variables
			const SphereData_Spectral &i_u,	///< prognostic variables
			const SphereData_Spectral &i_v,	///< prognostic variables

			SphereData_Spectral &o_h_t,	///< time updates
			SphereData_Spectral &o_u_t,	///< time updates
			SphereData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	void euler_timestep_update_coriolis_and_nonlinear(
			SphereData_Spectral &io_phi,		///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


private:
	void euler_timestep_update(
			const SphereData_Spectral &i_phi,	///< prognostic variables
			const SphereData_Spectral &i_vort,	///< prognostic variables
			const SphereData_Spectral &i_div,	///< prognostic variables

			SphereData_Spectral &o_phi_t,	///< time updates
			SphereData_Spectral &o_vort_t,	///< time updates
			SphereData_Spectral &o_div_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	SWE_Sphere_TS_lg_erk_lc_n_t_erk(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			SphereData_Spectral &io_phi,	///< prognostic variables
			SphereData_Spectral &io_vort,	///< prognostic variables
			SphereData_Spectral &io_div,	///< prognostic variables

			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_lg_erk_lc_n_t_erk();

        const SphereData_Spectral& getTopography() const
         {
	   return phi_topo;
	 }
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
