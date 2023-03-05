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
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "PDESWEPlaneTS_BaseInterface.hpp"
#include "SWE_Plane_TS_l_rexi.hpp"


class SWE_Plane_TS_l_rexi_na_sl_nr_etdrk	: public PDESWEPlaneTS_BaseInterface
{
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

	sweet::SemiLagrangianPlaneData semiLagrangian;
	sweet::PlaneDataSampler sampler2D;

	//Previous values (t_n-1)
	sweet::PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	sweet::ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	sweet::ScalarDataArray posx_d, posy_d;

	int timestepping_order;
	bool use_only_linear_divergence;


public:
	bool setup(
			sweet::PlaneOperators *io_ops
	);

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

	void euler_timestep_update_nonlinear(
			const sweet::PlaneData_Spectral &i_h,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_u,	///< prognostic variables
			const sweet::PlaneData_Spectral &i_v,	///< prognostic variables

			sweet::PlaneData_Spectral &o_h_t,	///< time updates
			sweet::PlaneData_Spectral &o_u_t,	///< time updates
			sweet::PlaneData_Spectral &o_v_t,	///< time updates

			double i_timestamp
	);


	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

#if ( SWEET_PARAREAL && SWEET_PARAREAL_PLANE ) || ( SWEET_XBRAID && SWEET_XBRAID_PLANE )
	void set_previous_solution(
				sweet::PlaneData_Spectral &i_h_prev,
				sweet::PlaneData_Spectral &i_u_prev,
				sweet::PlaneData_Spectral &i_v_prev
	) override
	{
		if (shackDict.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		h_prev = i_h_prev;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}
#endif

	virtual ~SWE_Plane_TS_l_rexi_na_sl_nr_etdrk() {}
};

#endif
