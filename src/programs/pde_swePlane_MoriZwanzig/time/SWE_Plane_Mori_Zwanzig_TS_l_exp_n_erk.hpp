/*
 * SWE_Plane_TS_l_exp_n_erk.hpp
 *
 *  Created on: 11 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_EXP_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_EXP_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "SWE_Plane_Mori_Zwanzig_TS_l_exp.hpp"
#include "SWE_Plane_Mori_Zwanzig_TS_n_erk.hpp"

#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "../PDESWEPlaneMoriZwanzig_Projection.hpp"

class SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk	: public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	int timestepping_order_nonlinear_P;
	int timestepping_order_nonlinear_Q;

	bool use_only_linear_divergence;

	sweet::TimesteppingExplicitRKPlaneData timestepping_rk_P;

	SWE_Plane_Mori_Zwanzig_TS_l_exp ts_l_exp;
	SWE_Plane_Mori_Zwanzig_TS_n_erk ts_n_erk;

	PDESWEPlaneMoriZwanzigProjection* projection = nullptr;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

public:
	bool setup(
			sweet::PlaneOperators *io_ops
	);

	void runTimestep(
			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SQ,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_FQ,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_FQ,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_exp_n_erk() {}
};

#endif
