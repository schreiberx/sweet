/*
 *  Created on: 30 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_PLANE_REXI_ADV_PLANE_TS_NA_SL_HPP_

#include <limits>
#include <sweet/ErrorBase.hpp>
#include <sweet/plane/Plane.hpp>

#include "PDEAdvPlaneTS_interface.hpp"


#include <sweet/shacks/ShackDictionary.hpp>
#include <sweet/shacksShared/ShackTimestepControl.hpp>
#include <sweet/shacksShared/ShackPlaneDataOps.hpp>
#include "ShackPDEAdvectionPlaneTimeDiscretization.hpp"
#include "../benchmarks/ShackPDEAdvectionPlaneBenchmarks.hpp"

#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>


class PDEAdvPlaneTS_na_sl	: public PDEAdvPlaneTS_interface
{
	sweet::ErrorBase error;

	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	ShackPDEAdvectionPlaneTimeDiscretization *shackTimeDisc;
	ShackPDEAdvectionPlaneBenchmarks *shackBenchmark;

	sweet::PlaneOperators &op;

	int timestepping_order;

	sweet::PlaneDataSampler sampler2D;
	sweet::PlaneDataSemiLagrangian semiLagrangian;

	sweet::PlaneData_Spectral prog_u_prev, prog_v_prev;


	/**
	 * Position in physical space given in longitude/latitude angles
	 */
	sweet::ScalarDataArray posx_a, posy_a;

public:
	PDEAdvPlaneTS_na_sl(
			sweet::ShackDictionary &io_shackDict,
			sweet::PlaneOperators &i_op
		);

private:
	void _setup(
			int i_order	///< order of RK time stepping method
	);

public:
	void run_timestep(
			sweet::PlaneData_Spectral &io_phi,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	virtual ~PDEAdvPlaneTS_na_sl();
};

#endif
