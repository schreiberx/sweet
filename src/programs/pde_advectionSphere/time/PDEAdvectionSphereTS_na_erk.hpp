/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_ERK_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_ERK_HPP_

#include "PDEAdvectionSphereTS_BaseInterface.hpp"
#include <limits>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingExplicitRKSphereData.hpp>
#include "../PDEAdvectionSphereBenchmarksCombined.hpp"



class PDEAdvectionSphereTS_na_erk	:
		public PDEAdvectionSphereTS_BaseInterface
{
	int timestepping_order;

	sweet::TimesteppingExplicitRKSphereData timestepping_rk;

public:
	bool testImplementsTimesteppingMethod(
			const std::string &i_timestepping_method
	);

	std::string getStringId();

	void printImplementedTimesteppingMethods(
			std::ostream &o_ostream,
			const std::string &i_prefix
	);


private:
	void euler_timestep_update(
			const sweet::SphereData_Spectral &i_prognostic_field,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			sweet::SphereData_Spectral &o_prognostic_field,	///< time updates

			double i_simulation_timestamp
	);

public:
	bool setup(
			sweet::SphereOperators *io_ops
	);

	void runTimestep(
			std::vector<sweet::SphereData_Spectral> &io_prognostic_fields,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	virtual ~PDEAdvectionSphereTS_na_erk();
};

#endif
