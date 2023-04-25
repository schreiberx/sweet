/*
 * SWE_Plane_TS_l_direct.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
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

#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>
#include "ShackPDESWEPlaneMoriZwanzigTimeDiscretization.hpp"
#include "../ShackPDESWEPlaneMoriZwanzig.hpp"
#include "../PDESWEPlaneMoriZwanzig_NormalModes.hpp"

///#include "../../pde_swePlane/time/PDESWEPlaneTS_BaseInterface.hpp"
#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"



class SWE_Plane_Mori_Zwanzig_TS_l_direct	:
		public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	typedef double T;

	sweet::ExpFunctions<T> expFunctions;

	///PDESWEPlaneMoriZwanzigNormalModes normal_modes;
	SWE_Plane_NormalModes normal_modes;

	sweet::PlaneDataGridMapping planeDataGridMapping;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

public:

	void runTimestep(
////			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,		///< prognostic variables



			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	void run_timestep_cgrid(
////			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_u,		///< prognostic variables
////			sweet::PlaneData_Spectral &io_v,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,		///< prognostic variables


			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid(
////			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,		///< prognostic variables


			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_agrid_planedata(
////			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
////			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,		///< prognostic variables


			double i_dt,
			double i_simulation_timestamp
	);

/////	void run_timestep_agrid_planedatacomplex(
/////			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
/////			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
/////			sweet::PlaneData_Spectral &io_v,	///< prognostic variables
/////
/////			double i_dt,
/////			double i_simulation_timestamp
/////	);

	bool setup(
			sweet::PlaneOperators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::PlaneOperators *io_ops
	);

	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_direct() {}
};

#endif
