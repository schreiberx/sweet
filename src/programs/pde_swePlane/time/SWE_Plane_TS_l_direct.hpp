/*
 * SWE_Plane_TS_l_direct.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_TS_L_DIRECT_HPP_

#include <sweet/expIntegration/ExpFunctions.hpp>
#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataGridMapping.hpp>
#include <sweet/core/plane/PlaneStaggering.hpp>

#include "SWE_Plane_TS_interface.hpp"

#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWEPlaneTimeDiscretization.hpp"
#include "../ShackPDESWEPlane.hpp"


class SWE_Plane_TS_l_direct	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackExpIntegration *shackExpIntegration;
	ShackPDESWEPlaneTimeDiscretization *shackTimeDisc;
	ShackPDESWEPlane *shackPDESWEPlane;

	sweet::PlaneOperators &op;


	int phi_id;

#if SWEET_QUADMATH && 0
	typedef __float128 T;
#else
	typedef double T;
#endif

	sweet::ExpFunctions<T> rexiFunctions;

	sweet::PlaneDataGridMapping planeDataGridMapping;

	// Precompute Z = Z(k1, k2, Dt) = Q * e^{\Lambda*Dt} * Q^{-1}
	std::vector<std::vector<std::array<std::array<std::complex<T>, 3>, 3>>> Z;  // Z[k1][k2][0,1,2][0,1,2];
	double dt_precompute_phin= 0.; // store dt for which Z has been precomputed in order to check if it is necessary to recompute it.


public:
	SWE_Plane_TS_l_direct(
			sweet::ShackDictionary *shackDict,
			sweet::PlaneOperators &i_op
		);


	void run_timestep(
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

	void setup(
			const std::string &i_function_name = "phi0"
	);

	virtual ~SWE_Plane_TS_l_direct();
};

#endif
