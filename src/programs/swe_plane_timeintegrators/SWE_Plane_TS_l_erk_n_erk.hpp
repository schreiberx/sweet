/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_ERK_N_ERK_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>

#include "SWE_Plane_TS_interface.hpp"



class SWE_Plane_TS_l_erk_n_erk	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int timestepping_order;
	int timestepping_order2;
	bool use_only_linear_divergence;

	PlaneDataTimesteppingExplicitRK timestepping_rk_linear;
	PlaneDataTimesteppingExplicitRK timestepping_rk_nonlinear;

private:
	void euler_timestep_update_linear(
			const PlaneData_Spectral &i_h,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_t,	///< time updates
			PlaneData_Spectral &o_u_t,	///< time updates
			PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


private:
	void euler_timestep_update_nonlinear(
			const PlaneData_Spectral &i_h,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_t,	///< time updates
			PlaneData_Spectral &o_u_t,	///< time updates
			PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	SWE_Plane_TS_l_erk_n_erk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2,
			bool i_use_only_linear_divergence = false
	);

	void run_timestep(
			PlaneData_Spectral &io_h_pert,	///< prognostic variables
			PlaneData_Spectral &io_u,		///< prognostic variables
			PlaneData_Spectral &io_v,		///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);





	virtual ~SWE_Plane_TS_l_erk_n_erk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_LN_ERK_HPP_ */
