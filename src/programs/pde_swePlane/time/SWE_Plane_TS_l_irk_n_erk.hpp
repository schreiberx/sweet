/*
 * SWE_Plane_TS_l_irk_n_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_IRK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_IRK_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "SWE_Plane_TS_interface.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"


class SWE_Plane_TS_l_irk_n_erk	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	sweet::PlaneOperators &op;

	int timestepping_order_linear;
	int timestepping_order_nonlinear;

	bool use_only_linear_divergence;

	sweet::PlaneDataTimesteppingExplicitRK timestepping_rk;
	SWE_Plane_TS_l_irk ts_l_irk;


public:
	void euler_timestep_update_nonlinear(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t	///< time updates
	);


public:
	SWE_Plane_TS_l_irk_n_erk(
			sweet::ShackDictionary *shackDict,
			sweet::PlaneOperators &i_op
		);

	void setup(
			int i_l_order,
			int i_n_order,

			bool i_use_only_linear_divergence
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

        SWE_Plane_TS_l_irk& get_implicit_timestepper() {return ts_l_irk;}

	virtual ~SWE_Plane_TS_l_irk_n_erk();
};

#endif
