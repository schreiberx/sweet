/*
 * SWE_Plane_TS_l_cn_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Modif on: 11 July 2017
 *  	Author: Pedro Peixoto <pedrosp@ime.usp.br>
 *
 *  Name: l_cn_na_sl_nd_settls
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_CN_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>
#include <sweet/core/time/TimesteppingSemiLagrangianPlaneData.hpp>


#include "PDESWEPlaneTS_BaseInterface.hpp"

class SWE_Plane_TS_l_cn_na_sl_nr_settls	: public PDESWEPlaneTS_BaseInterface
{
	bool use_only_linear_divergence;

	sweet::TimesteppingSemiLagrangianPlaneData semiLagrangian;
	sweet::PlaneDataSampler sampler2D;

	sweet::PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	sweet::ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	sweet::ScalarDataArray posx_d, posy_d;

public:
	bool setup(
			sweet::PlaneOperators *io_ops
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);



	void helmholtz_spectral_solver(
			double i_kappa,
			double i_gh0,
			const sweet::PlaneData_Spectral &i_rhs,
			sweet::PlaneData_Spectral &io_x
	)
	{
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		sweet::PlaneData_Spectral laplacian = -i_gh0*ops->diff2_c_x -i_gh0*ops->diff2_c_y;
		sweet::PlaneData_Spectral lhs = laplacian.spectral_addScalarAll(i_kappa);

		io_x = i_rhs.spectral_div_element_wise(lhs);
#else
		SWEETError("Cannot use helmholtz_spectral_solver if spectral space not enable in compilation time");
#endif
	}

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

	virtual ~SWE_Plane_TS_l_cn_na_sl_nr_settls();
};

#endif
