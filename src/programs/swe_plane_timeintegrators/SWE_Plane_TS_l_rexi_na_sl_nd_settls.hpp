/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "../swe_plane_timeintegrators/SWE_Plane_TS_interface.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_direct.hpp"
#include "../swe_plane_timeintegrators/SWE_Plane_TS_l_rexi.hpp"



class SWE_Plane_TS_l_rexi_na_sl_nd_settls	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	SWE_Plane_TS_l_rexi ts_l_rexi;

	int with_linear_div_only;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	//Previous values (t_n-1)
	PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	bool use_only_linear_divergence;

public:
	SWE_Plane_TS_l_rexi_na_sl_nd_settls(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			//REXI_SimulationVariables &i_rexi,
			//int i_with_nonlinear
			bool use_only_linear_divergence
	);

	void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);

#if ( SWEET_PARAREAL && SWEET_PARAREAL_PLANE ) || ( SWEET_XBRAID && SWEET_XBRAID_PLANE )
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
#endif



	virtual ~SWE_Plane_TS_l_rexi_na_sl_nd_settls();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_ */
