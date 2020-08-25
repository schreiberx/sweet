/*
 * SWE_Plane_TS_PFASST_ln_erk.hpp
 *
 *  Created on: 29 May 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_SWE_SPHERE_TS_PFASST_INTERFACE_PFASST_HPP_
#define SRC_PROGRAMS_SWE_SPHERE_TS_PFASST_INTERFACE_PFASST_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>


class SWE_Sphere_TS_PFASST_interface
{
public:
	virtual void run_timestep_nonpert(
			SphereData_Spectral &io_phi_pert,	///< prognostic variables
			SphereData_Spectral &io_vrt,		///< prognostic variables
			SphereData_Spectral &io_div,		///< prognostic variables

			double i_fixed_dt,		
			double i_simulation_timestamp
	) = 0;
};

#endif /* SRC_PROGRAMS_SWE_PLANE_REXI_SWE_PLANE_TS_PFASST_LN_ERK_HPP_ */
