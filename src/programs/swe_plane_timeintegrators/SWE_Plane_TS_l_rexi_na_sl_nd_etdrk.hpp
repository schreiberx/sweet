/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_etdrk.hpp
 *
 *  Created on: 09 Oct 2017
 *      Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Changelog:
 *      based on Martin Schreiber ETD timestepper
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_ETDRK_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_ETDRK_HPP_

#include <limits>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "../swe_plane_timeintegrators/SWE_Plane_TS_interface.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_rexi.hpp"


class SWE_Plane_TS_l_rexi_na_sl_nd_etdrk	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	SWE_Plane_TS_l_rexi ts_phi0_rexi;
	SWE_Plane_TS_l_rexi ts_phi1_rexi;
	SWE_Plane_TS_l_rexi ts_phi2_rexi;

	SWE_Plane_TS_l_rexi ts_ups0_rexi;
	SWE_Plane_TS_l_rexi ts_ups1_rexi;
	SWE_Plane_TS_l_rexi ts_ups2_rexi;
	SWE_Plane_TS_l_rexi ts_ups3_rexi;

	SWE_Plane_TS_l_rexi ts_psi1_rexi;
	SWE_Plane_TS_l_rexi ts_psi2_rexi;
	SWE_Plane_TS_l_rexi ts_psi3_rexi;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	//Previous values (t_n-1)
	PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	int timestepping_order;
	bool use_only_linear_divergence;


public:
	SWE_Plane_TS_l_rexi_na_sl_nd_etdrk(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);


	void setup(
		int i_timestepping_order,
		bool i_use_only_linear_divergence
	);


	void euler_timestep_update_nonlinear(
			const PlaneData_Spectral &i_h,	///< prognostic variables
			const PlaneData_Spectral &i_u,	///< prognostic variables
			const PlaneData_Spectral &i_v,	///< prognostic variables

			PlaneData_Spectral &o_h_t,	///< time updates
			PlaneData_Spectral &o_u_t,	///< time updates
			PlaneData_Spectral &o_v_t,	///< time updates

			double i_timestamp
	);


	void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

	void set_previous_solution(
				PlaneData_Spectral &i_h_prev,
				PlaneData_Spectral &i_u_prev,
				PlaneData_Spectral &i_v_prev
	) override
	{
		if (simVars.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		h_prev = i_h_prev;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}



	virtual ~SWE_Plane_TS_l_rexi_na_sl_nd_etdrk();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_l_rexi_na_sl_nd_etdrk_HPP_ */
