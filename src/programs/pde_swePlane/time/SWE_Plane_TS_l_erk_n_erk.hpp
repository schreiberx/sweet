/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_ERK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_ERK_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>

#include "SWE_Plane_TS_interface.hpp"



class SWE_Plane_TS_l_erk_n_erk	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	sweet::PlaneOperators &op;

	int timestepping_order;
	int timestepping_order2;
	bool use_only_linear_divergence;

	sweet::PlaneDataTimesteppingExplicitRK timestepping_rk_linear;
	sweet::PlaneDataTimesteppingExplicitRK timestepping_rk_nonlinear;

private:
	void euler_timestep_update_linear(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


private:
	void euler_timestep_update_nonlinear(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp
	);


public:
	SWE_Plane_TS_l_erk_n_erk(
			sweet::ShackDictionary *shackDict,
			sweet::PlaneOperators &i_op
		);

	void setup(
			int i_order,	///< order of RK time stepping method
			int i_order2,
			bool i_use_only_linear_divergence = false
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v,		///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);





	virtual ~SWE_Plane_TS_l_erk_n_erk();
};

#endif
