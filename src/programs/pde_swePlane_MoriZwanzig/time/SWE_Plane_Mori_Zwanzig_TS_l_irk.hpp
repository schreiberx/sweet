/*
 * SWE_Plane_Mori_Zwanzig_TS_l_irk.hpp
 *
 *  Created on: 18 May 2023
 *  Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_XX_SWE_PLANE_MORI_ZWANZIG_SWE_PLANE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_XX_SWE_PLANE_MORI_ZWANZIG_SWE_PLANE_TS_L_IRK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>

#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "SWE_Plane_Mori_Zwanzig_TS_l_direct.hpp"

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#endif


class SWE_Plane_Mori_Zwanzig_TS_l_irk	: public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	sweet::PlaneOperatorsComplex opComplex;
#endif

	int timestepping_order_P;
	int timestepping_order_Q;

public:
	bool setup(
			sweet::PlaneOperators *io_ops,
			int i_order_P,
			int i_order_Q
	);

	bool setup(
			sweet::PlaneOperators *io_ops
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



	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_irk() {}
};

#endif
