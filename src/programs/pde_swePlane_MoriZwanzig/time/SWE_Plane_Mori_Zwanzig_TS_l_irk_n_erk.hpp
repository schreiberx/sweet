/*
 * SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk.hpp
 *
 *  Created on: 18 May 2023
 *  Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_IRK_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_MORI_ZWANZIG_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_IRK_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "SWE_Plane_Mori_Zwanzig_TS_l_irk.hpp"
#include "SWE_Plane_Mori_Zwanzig_TS_n_erk.hpp"


class SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk	:
		public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	int timestepping_order_linear_P;
	int timestepping_order_nonlinear_P;
	int timestepping_order_linear_Q;
	int timestepping_order_nonlinear_Q;

	bool use_only_linear_divergence;

	SWE_Plane_Mori_Zwanzig_TS_l_irk ts_l_irk;
	SWE_Plane_Mori_Zwanzig_TS_n_erk ts_n_erk;

	std::string equation;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);


public:
	bool setup(
			sweet::PlaneOperators *io_ops,
			std::string i_equation
	);

	void runTimestep(
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

	SWE_Plane_Mori_Zwanzig_TS_l_irk& get_implicit_timestepper() {return ts_l_irk;}

	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_irk_n_erk() {}
};

#endif
