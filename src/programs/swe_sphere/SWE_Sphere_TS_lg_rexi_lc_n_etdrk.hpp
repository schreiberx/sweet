/*
 * SWE_Sphere_TS_lg_rexi_lf_n_etdrk.hpp
 *
 *  Created on: 10 Sept. 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_Sphere_TS_lg_rexi_lf_n_etdrk_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_Sphere_TS_lg_rexi_lf_n_etdrk_HPP_

#include <limits>
#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataPhysicalTimesteppingExplicitRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include "SWE_Sphere_TS_interface.hpp"

#include "SWE_Sphere_TS_l_rexi.hpp"
#include "SWE_Sphere_TS_lg_erk_lc_n_erk.hpp"


class SWE_Sphere_TS_lg_rexi_lc_n_etdrk	: public SWE_Sphere_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators &op;

	SWE_Sphere_TS_lg_erk_lc_n_erk ts_lg_erk_lf_n_erk;

	SWE_Sphere_TS_l_rexi ts_phi0_rexi;
	SWE_Sphere_TS_l_rexi ts_phi1_rexi;
	SWE_Sphere_TS_l_rexi ts_phi2_rexi;

	SWE_Sphere_TS_l_rexi ts_ups0_rexi;
	SWE_Sphere_TS_l_rexi ts_ups1_rexi;
	SWE_Sphere_TS_l_rexi ts_ups2_rexi;
	SWE_Sphere_TS_l_rexi ts_ups3_rexi;

	int timestepping_order;
	int timestepping_order2;

private:
	void euler_timestep_update_linear(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double i_max_timestamp
	);



private:
	void euler_timestep_update_nonlinear(
			const SphereData &i_h,	///< prognostic variables
			const SphereData &i_u,	///< prognostic variables
			const SphereData &i_v,	///< prognostic variables

			SphereData &o_h_t,	///< time updates
			SphereData &o_u_t,	///< time updates
			SphereData &o_v_t,	///< time updates

			double i_max_timestamp
	);


public:
	SWE_Sphere_TS_lg_rexi_lc_n_etdrk(
			SimulationVariables &i_simVars,
			SphereOperators &i_op
		);

	void setup(
			REXI_SimulationVariables &i_rexi,
			int i_timestepping_order,
			int i_timestepping_order2,
			double i_timestep_size
	);

	void run_timestep(
			SphereData &io_h,	///< prognostic variables
			SphereData &io_u,	///< prognostic variables
			SphereData &io_v,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	virtual ~SWE_Sphere_TS_lg_rexi_lc_n_etdrk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
