/*
 * SWE_Plane_Mori_Zwanzig_TS_n_erk.hpp
 *
 *  Created on: 11 Apr 2023
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_N_ERK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_TIMEINTEGRATORS_SWE_PLANE_MORI_ZWANZIG_TS_L_N_ERK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include "../../pde_swePlane/time/SWE_Plane_TS_l_rexi.hpp"

#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"
#include "../PDESWEPlaneMoriZwanzig_Projection.hpp"

class SWE_Plane_Mori_Zwanzig_TS_n_erk	: public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	int timestepping_order_nonlinear_P;
	int timestepping_order_nonlinear_Q;

	// runge kutta data storages
	PlaneData_Spectral** RK_h_SQ_t = nullptr;
	PlaneData_Spectral** RK_u_SQ_t = nullptr;
	PlaneData_Spectral** RK_v_SQ_t = nullptr;
	PlaneData_Spectral** RK_h_FQ_t = nullptr;
	PlaneData_Spectral** RK_u_FQ_t = nullptr;
	PlaneData_Spectral** RK_v_FQ_t = nullptr;


	bool use_only_linear_divergence;

	sweet::TimesteppingExplicitRKPlaneData timestepping_rk_P;

	PDESWEPlaneMoriZwanzigProjection* projection = nullptr;

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

private:
	void euler_timestep_update_nonlinear_P(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_timestamp
	);

private:
	void euler_timestep_update_nonlinear_Q(
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


	void SWE_Plane_Mori_Zwanzig_TS_n_erk:setupBuffers(
			const PlaneData_Config *i_planeDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
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

	void runTimestep_P(
			sweet::PlaneData_Spectral &io_h_pert_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_SP,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_SP,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	void runTimestep_Q(
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
