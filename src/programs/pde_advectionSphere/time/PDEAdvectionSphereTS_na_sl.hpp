/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_


#include "../PDEAdvectionSphereBenchmarksCombined.hpp"
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <limits>
#include <sweet/core/shacks/ShackDictionary.hpp>

#include <sweet/core/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/time/TimesteppingSemiLagrangianSphereData.hpp>
#include "../time/PDEAdvectionSphereTS_BaseInterface.hpp"


class PDEAdvectionSphereTS_na_sl	:
		public PDEAdvectionSphereTS_BaseInterface
{
	int timestepping_order;

	sweet::TimesteppingSemiLagrangianSphereData semiLagrangian;
	sweet::SphereOperators_Sampler_SphereDataPhysical sphereSampler;

	sweet::SphereData_Physical U_u_prev, U_v_prev;

public:
	bool testImplementsTimesteppingMethod(const std::string &i_timestepping_method);

	std::string getStringId();

	void printImplementedTimesteppingMethods(
			std::ostream &o_ostream,
			const std::string &i_prefix
	);


public:
	bool setup(
			sweet::SphereOperators *io_ops
	);

private:
	void interpolate_departure_point_vec_3d(
			const sweet::SphereData_Spectral &i_u,
			const sweet::SphereData_Spectral &i_v,
			const sweet::SphereData_Spectral &i_w,

			const sweet::ScalarDataArray &i_pos_lon_D,
			const sweet::ScalarDataArray &i_pos_lat_D,

			sweet::SphereData_Spectral &o_u,
			sweet::SphereData_Spectral &o_v,
			sweet::SphereData_Spectral &o_w
	);

	void interpolate_departure_point_vec_uv(
			const sweet::SphereData_Physical &i_u,
			const sweet::SphereData_Physical &i_v,

			const sweet::ScalarDataArray &i_pos_lon_D,
			const sweet::ScalarDataArray &i_pos_lat_D,

			sweet::SphereData_Physical &o_u,
			sweet::SphereData_Physical &o_v
	);

	void runTimestep(
			std::vector<sweet::SphereData_Spectral> &io_prog_fields,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);


	void run_timestep_1(
			sweet::SphereData_Spectral &io_prognostic_field,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);


	void run_timestep_2(
			std::vector<sweet::SphereData_Spectral> &io_prog_fields,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	void run_timestep_3(
			std::vector<sweet::SphereData_Spectral> &io_prog_fields,	///< prognostic variables
			sweet::SphereData_Physical &io_u,
			sweet::SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp
	);



	virtual ~PDEAdvectionSphereTS_na_sl();
};

#endif
