/*
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_

#include "PDEAdvectionPlaneTS_BaseInterface.hpp"
#include <limits>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/plane/Plane.hpp>

#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/time/TimesteppingSemiLagrangianPlaneData.hpp>
#include "../benchmarks/ShackPDEAdvectionPlaneBenchmarks.hpp"


class PDEAdvectionPlaneTS_na_sl	:
		public PDEAdvectionPlaneTS_BaseInterface
{
	int timestepping_order;

	sweet::PlaneDataSampler sampler2D;
	sweet::TimesteppingSemiLagrangianPlaneData semiLagrangian;

	sweet::PlaneData_Spectral prog_u_prev, prog_v_prev;


	/**
	 * Position in physical space given in longitude/latitude angles
	 */
	sweet::ScalarDataArray posx_a, posy_a;

public:
	PDEAdvectionPlaneTS_na_sl();

	~PDEAdvectionPlaneTS_na_sl();

	bool setup(sweet::PlaneOperators *io_ops);

private:
	void _setup();

public:
	void runTimestep(
			sweet::PlaneData_Spectral &io_phi,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);
};

#endif
