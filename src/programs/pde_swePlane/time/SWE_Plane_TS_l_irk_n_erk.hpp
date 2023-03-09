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
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "PDESWEPlaneTS_BaseInterface.hpp"
#include "SWE_Plane_TS_l_irk.hpp"
#include "SWE_Plane_TS_l_irk.hpp"


class SWE_Plane_TS_l_irk_n_erk	:
		public PDESWEPlaneTS_BaseInterface
{
	int timestepping_order_linear;
	int timestepping_order_nonlinear;

	bool use_only_linear_divergence;

	sweet::TimesteppingExplicitRKPlaneData timestepping_rk;
	SWE_Plane_TS_l_irk ts_l_irk;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

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
	bool setup(
			sweet::PlaneOperators *io_ops
	);

	void runTimestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	SWE_Plane_TS_l_irk& get_implicit_timestepper() {return ts_l_irk;}

	virtual ~SWE_Plane_TS_l_irk_n_erk() {}
};

#endif
