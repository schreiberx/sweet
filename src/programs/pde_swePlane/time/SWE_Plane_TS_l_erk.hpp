/*
 * SWE_Plane_TS_l_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include "PDESWEPlaneTS_BaseInterface.hpp"



class SWE_Plane_TS_l_erk	:
		public PDESWEPlaneTS_BaseInterface
{
	int timestepping_order;

	// Sampler
	sweet::PlaneDataTimesteppingExplicitRK timestepping_rk;

public:
	void euler_timestep_update(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_simulation_timestamp = -1
	);

public:
	bool setup(
			sweet::PlaneOperators *io_ops,
			int i_order
		);


	bool setup(
			sweet::PlaneOperators *io_ops
		);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Plane_TS_l_erk() {}
};

#endif
