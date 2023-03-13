/*
 * SWE_Plane_TS_l_rexi_na_sl_nd_settls.hpp
 *
 *  Created on: 29 May 2017
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Changelog:
 *  	2017-05-29: Based on source swe_plane.cpp
 *					which was also written by Pedro Peixoto
 */

#ifndef SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_
#define SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_L_REXI_NA_SL_ND_SETTLS_HPP_

#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>
#include <sweet/core/plane/PlaneDataSampler.hpp>
#include <sweet/core/time/TimesteppingSemiLagrangianPlaneData.hpp>
#include <sweet/core/time/TimesteppingExplicitRKPlaneData.hpp>

#include "PDESWEPlaneTS_BaseInterface.hpp"
#include "SWE_Plane_TS_l_direct.hpp"
#include "SWE_Plane_TS_l_rexi.hpp"



class SWE_Plane_TS_l_rexi_na_sl_nr_settls	: public PDESWEPlaneTS_BaseInterface
{



	SWE_Plane_TS_l_rexi ts_l_rexi;

	//int with_linear_div_only;

	sweet::TimesteppingSemiLagrangianPlaneData semiLagrangian;
	sweet::PlaneDataSampler sampler2D;

	//Previous values (t_n-1)
	sweet::PlaneData_Spectral h_prev, u_prev, v_prev;

	// Arrival points for semi-lag
	sweet::ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	sweet::ScalarDataArray posx_d, posy_d;

	bool use_only_linear_divergence;

public:
	bool setup(
			sweet::PlaneOperators *io_ops
	);

	bool shackRegistration(
			sweet::ShackDictionary *io_shackDict
	);

	void runTimestep(
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

	virtual ~SWE_Plane_TS_l_rexi_na_sl_nr_settls() {}
};

#endif
