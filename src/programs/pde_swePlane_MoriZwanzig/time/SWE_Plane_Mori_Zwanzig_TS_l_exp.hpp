/*
 * SWE_Plane_TS_ln_erk.hpp
 *
 *  Created on: 14 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TS_L_EXP_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TS_L_EXP_HPP_

#include <sweet/core/defaultPrecompilerValues.hpp>
#include <limits>
#include <string>
#include <complex>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneData_SpectralComplex.hpp>
#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/expIntegration/ShackExpIntegration.hpp>

///#include "../../pde_swePlane/time/PDESWEPlaneTS_BaseInterface.hpp"
#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "SWE_Plane_Mori_Zwanzig_TS_l_direct.hpp"

#if SWEET_MPI
#	include <mpi.h>
#endif

class SWE_Plane_Mori_Zwanzig_TS_l_exp	:
		public PDESWEPlaneMoriZwanzigTS_BaseInterface
{

public:
	/// final time step
	bool final_timestep;

	/// use direct solution instead of REXI
	bool exp_use_direct_solution;

	/// Direct solution for linear parts
	SWE_Plane_Mori_Zwanzig_TS_l_direct ts_l_direct;

public:
	bool setup(
			sweet::PlaneOperators *io_ops,
			const std::string &i_function_name
	);

	bool setup(
			sweet::PlaneOperators *io_ops
	);

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

	void runTimestep(
///			const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
///			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
///			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

///			sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
///			sweet::PlaneData_Spectral &o_u,	///< prognostic variables
///			sweet::PlaneData_Spectral &o_v,	///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_SP,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_SP,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_SP,		///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_SQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_SQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_SQ,		///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_FQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_FQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_FQ,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_FQ,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);


	void run_timestep_real(
//			const sweet::PlaneData_Spectral &i_h_pert,	///< prognostic variables
//			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
//			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables
//
//			sweet::PlaneData_Spectral &o_h_pert,	///< prognostic variables
//			sweet::PlaneData_Spectral &o_u,	///< prognostic variables
//			sweet::PlaneData_Spectral &o_v,	///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_SP,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_SP,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_SP,		///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_SQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_SQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_SQ,		///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_pert_FQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_FQ,		///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_FQ,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_SP,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_SP,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_SQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_SQ,		///< prognostic variables

			sweet::PlaneData_Spectral &o_h_pert_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_u_FQ,		///< prognostic variables
			sweet::PlaneData_Spectral &o_v_FQ,		///< prognostic variables

			double i_dt,
			double i_simulation_timestamp
	);

	void runTimestep(
///			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
///			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
///			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

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



	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_exp();
};

#endif
