/*
 * SWE_Plane_TS_l_cn_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
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
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include "SWE_Plane_TS_interface.hpp"



class SWE_Plane_TS_l_cn_na_sl_nd_settls	: public SWE_Plane_TS_interface
{
	SimulationVariables &simVars;
	PlaneOperators &op;

	int with_nonlinear;

	PlaneDataSemiLagrangian semiLagrangian;
	PlaneDataSampler sampler2D;

	PlaneData h_prev, u_prev, v_prev;

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
			int i_with_nonlinear
	);

	void run_timestep(
			PlaneData &io_h,	///< prognostic variables
			PlaneData &io_u,	///< prognostic variables
			PlaneData &io_v,	///< prognostic variables

			double i_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	);



	void helmholtz_spectral_solver(
			double i_kappa,
			double i_gh0,
			const PlaneData &i_rhs,
			PlaneData &io_x
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PlaneData laplacian = -i_gh0*op.diff2_c_x -i_gh0*op.diff2_c_y;
		PlaneData lhs = laplacian.spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
#else
		FatalError("Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time");
#endif
	}


	virtual ~SWE_Plane_TS_l_cn_na_sl_nd_settls();
};

#endif /* SRC_PROGRAMS_SWE_Plane_TS_l_cn_na_sl_nd_settls_HPP_ */
