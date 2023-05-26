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

#include "PDESWEPlaneMoriZwanzigTS_BaseInterface.hpp"

class SWE_Plane_Mori_Zwanzig_TS_n_erk	: public PDESWEPlaneMoriZwanzigTS_BaseInterface
{
	int timestepping_order_nonlinear_P;
	int timestepping_order_nonlinear_Q;

	// runge kutta data storages
	sweet::PlaneData_Spectral** RK_h_A_t = nullptr;
	sweet::PlaneData_Spectral** RK_u_A_t = nullptr;
	sweet::PlaneData_Spectral** RK_v_A_t = nullptr;
	sweet::PlaneData_Spectral** RK_h_B_t = nullptr;
	sweet::PlaneData_Spectral** RK_u_B_t = nullptr;
	sweet::PlaneData_Spectral** RK_v_B_t = nullptr;


	bool use_only_linear_divergence;

	sweet::TimesteppingExplicitRKPlaneData timestepping_rk_P;

	std::string equation; // SF or PQ

public:
	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

private:
	void euler_timestep_update_bilinear(
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

			sweet::PlaneData_Spectral &o_h_A_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_A_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_A_t,	///< time updates

			sweet::PlaneData_Spectral &o_h_B_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_B_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_B_t,	///< time updates

			double i_timestamp
	);

private:
	void euler_timestep_update_nonlinear_SF(
			const sweet::PlaneData_Spectral &i_h_A,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_A,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_A,	///< prognostic variables

			const sweet::PlaneData_Spectral &i_h_B,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u_B,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v_B,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_A_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_A_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_A_t,	///< time updates

			sweet::PlaneData_Spectral &o_h_B_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_B_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_B_t,	///< time updates

			double i_timestamp
	);


	void setupBuffers(
			const sweet::PlaneData_Config *i_planeDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
	);

public:
	bool setup(
			sweet::PlaneOperators *io_ops,
			const sweet::PlaneData_Config *i_planeDataConfig,
			std::string i_equation		// SF, P or Q
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

	void runTimestep_SF(
			sweet::PlaneData_Spectral &io_h_pert_S,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_S,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_S,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_F,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_F,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_F,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	void runTimestep_full(
			sweet::PlaneData_Spectral &io_h_pert,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


private:
	void runTimestep_Q_SF(
			void (SWE_Plane_Mori_Zwanzig_TS_n_erk::*i_compute_euler_timestep_update)(
					const sweet::PlaneData_Spectral &i_h_A,
					const sweet::PlaneData_Spectral &i_u_A,
					const sweet::PlaneData_Spectral &i_v_A,

					const sweet::PlaneData_Spectral &i_h_B,
					const sweet::PlaneData_Spectral &i_u_B,
					const sweet::PlaneData_Spectral &i_v_B,

					sweet::PlaneData_Spectral &o_h_A,
					sweet::PlaneData_Spectral &o_u_A,
					sweet::PlaneData_Spectral &o_v_A,

					sweet::PlaneData_Spectral &o_h_B,
					sweet::PlaneData_Spectral &o_u_B,
					sweet::PlaneData_Spectral &o_v_B,

					double i_dt
			),

			sweet::PlaneData_Spectral &io_h_pert_A,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_A,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_A,	///< prognostic variables

			sweet::PlaneData_Spectral &io_h_pert_B,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u_B,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v_B,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);




public:
	virtual ~SWE_Plane_Mori_Zwanzig_TS_n_erk();
};

#endif
