/*
 * SWE_Plane_TS_l_rexi_n_erk.hpp
 *
 *  Created on: 11 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_REXI_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_REXI_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/time/TimesteppingExplicitRKBilinearPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "../../pde_swePlane/time/SWE_Plane_TS_l_rexi.hpp"

#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "../PDESWEPlaneMoriZwanzig_Projection.hpp"

class SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk	: public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	int timestepping_order_nonlinear_SP;
	int timestepping_order_nonlinear_SQ;
	int timestepping_order_nonlinear_FQ;

	bool use_only_linear_divergence;

	sweet::TimesteppingExplicitRKBilinearPlaneData timestepping_rk_SP;
	sweet::TimesteppingExplicitRKBilinearPlaneData timestepping_rk_SQ;
	sweet::TimesteppingExplicitRKBilinearPlaneData timestepping_rk_FQ;

	SWE_Plane_TS_l_rexi ts_l_rexi;

	PDESWEPlaneMoriZwanzigProjection* projection = nullptr;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

private:
	void euler_timestep_update_nonlinear(
			const sweet::PlaneData_Spectral &i_h_A,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_A,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_A,	///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_B,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_B,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_B,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_timestamp
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

	virtual ~SWE_Plane_Mori_Zwanzig_TS_l_rexi_n_erk() {}
};

#endif
