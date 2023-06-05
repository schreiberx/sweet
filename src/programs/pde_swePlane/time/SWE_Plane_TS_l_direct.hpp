/*
 * SWE_Plane_TS_l_direct.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataGridMapping.hpp>
#include <sweet/core/plane/PlaneStaggering.hpp>

#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/expIntegration/ExpFunction.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWEPlaneTimeDiscretization.hpp"
#include "../ShackPDESWEPlane.hpp"

#include "PDESWEPlaneTS_BaseInterface.hpp"


class SWE_Plane_TS_l_direct	:
		public PDESWEPlaneTS_BaseInterface
{
	typedef double T;

	sweet::ExpFunction<T> expFunctions;

	sweet::PlaneDataGridMapping planeDataGridMapping;

	// Precompute Z = Z(k1, k2, Dt) = Q * e^{\Lambda*Dt} * Q^{-1}
	std::vector<std::vector<std::array<std::array<std::complex<T>, 3>, 3>>> Z;  // Z[k1][k2][0,1,2][0,1,2];
	double dt_precompute_phin = 0.; // store dt for which Z has been precomputed in order to check if it is necessary to recompute it.


public:

	void runTimestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_cgrid(
			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid_planedata(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);



	void run_timestep_agrid_planedatacomplex(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	bool setup(
			sweet::PlaneOperators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::PlaneOperators *io_ops
	);

	virtual ~SWE_Plane_TS_l_direct() {}
};

#endif
