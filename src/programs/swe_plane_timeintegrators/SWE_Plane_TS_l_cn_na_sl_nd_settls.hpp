/*
 * SWE_Plane_TS_l_cn_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Modif on: 11 July 2017
 *  	Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Name: l_cn_na_sl_nd_settls
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataTimesteppingExplicitRK.hpp>

#include "../swe_plane_timeintegrators/SWE_Plane_TS_interface.hpp"



class SWE_Plane_TS_l_cn_na_sl_nd_settls	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	bool use_only_linear_divergence;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

public:
	SWE_Plane_TS_l_cn_na_sl_nd_settls(
			SimulationVariables &i_simVars,
			PlaneOperators &i_op
		);

	void setup(
			bool i_use_only_linear_divergence
	);

	void run_timestep(
			PlaneData_Spectral &io_h,	///< prognostic variables
			PlaneData_Spectral &io_u,	///< prognostic variables
			PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	void helmholtz_spectral_solver(
			double i_kappa,
			double i_gh0,
			const PlaneData_Spectral &i_rhs,
			PlaneData_Spectral &io_x
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData_Spectral laplacian = -i_gh0*op.diff2_c_x -i_gh0*op.diff2_c_y;
		PlaneData_Spectral lhs = laplacian.spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
#else
		SWEETError("Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time");
#endif
	}

	void set_previous_solution(
				PlaneData &i_h_prev,
				PlaneData &i_u_prev,
				PlaneData &i_v_prev
	) override
	{
		if (simVars.misc.verbosity > 5)
			std::cout << "set_previous_solution()" << std::endl;
		h_prev = i_h_prev;
		u_prev = i_u_prev;
		v_prev = i_v_prev;
	}


	virtual ~SWE_Plane_TS_l_cn_na_sl_nd_settls();
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_NA_SL_ND_SETTLS_HPP_ */
