/*
 * SWE_Plane_TS_l_irk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_XX_SWE_PLANE_REXI_SWE_PLANE_TS_L_IRK_HPP_
#define SRC_PROGRAMS_XX_SWE_PLANE_REXI_SWE_PLANE_TS_L_IRK_HPP_

#include <limits>
#include <sweet/core/plane/PlaneData_Spectral.hpp>
#include <sweet/core/plane/PlaneDataTimesteppingExplicitRK.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include <sweet/core/plane/PlaneOperators.hpp>

#include "SWE_Plane_TS_interface.hpp"

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#include <sweet/core/plane/PlaneOperatorsComplex.hpp>
#endif


class SWE_Plane_TS_l_irk	: public SWE_Plane_TS_interface
{
	sweet::ShackDictionary *shackDict;
	sweet::PlaneOperators &op;

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	sweet::PlaneOperatorsComplex opComplex;
#endif

	int timestepping_order;

public:
	SWE_Plane_TS_l_irk(
			sweet::ShackDictionary *io_shackDict,
			sweet::PlaneOperators &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);

	void run_timestep(
			sweet::PlaneData_Spectral &io_h,	///< prognostic variables
			sweet::PlaneData_Spectral &io_u,	///< prognostic variables
			sweet::PlaneData_Spectral &io_v,	///< prognostic variables

			double i_dt = 0,
			double i_simulation_timestamp = -1
	);


	virtual ~SWE_Plane_TS_l_irk();
};

#endif
