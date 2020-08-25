/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#ifndef SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_
#define SRC_PROGRAMS_ADV_SPHERE_REXI_ADV_SPHERE_TS_NA_SL_HPP_


#include "../advection_sphere_benchmarks/BenchmarksSphereAdvection.hpp"
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <limits>
#include <sweet/SimulationVariables.hpp>

#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/sphere/SphereTimestepping_SemiLagrangian.hpp>

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_interface.hpp"


class SphereAdvection_TS_na_sl	: public SphereAdvection_TS_interface
{
	SimulationVariables &simVars;
	SphereOperators_SphereData &op;

	int timestepping_order;

	SphereTimestepping_SemiLagrangian semiLagrangian;
	SphereOperators_Sampler_SphereDataPhysical &sphereSampler;

	SphereData_Physical U_u_prev, U_v_prev;

public:
	bool implements_timestepping_method(const std::string &i_timestepping_method);

	std::string string_id();

	void setup_auto();

	std::string get_help();

public:
	SphereAdvection_TS_na_sl(
			SimulationVariables &i_simVars,
			SphereOperators_SphereData &i_op
		);

	void setup(
			int i_order	///< order of RK time stepping method
	);


	void interpolate_departure_point_vec_3d(
			const SphereData_Spectral &i_u,
			const SphereData_Spectral &i_v,
			const SphereData_Spectral &i_w,

			const ScalarDataArray &i_pos_lon_D,
			const ScalarDataArray &i_pos_lat_D,

			SphereData_Spectral &o_u,
			SphereData_Spectral &o_v,
			SphereData_Spectral &o_w
	);

	void interpolate_departure_point_vec_uv(
			const SphereData_Physical &i_u,
			const SphereData_Physical &i_v,

			const ScalarDataArray &i_pos_lon_D,
			const ScalarDataArray &i_pos_lat_D,

			SphereData_Physical &o_u,
			SphereData_Physical &o_v
	);

	void run_timestep(
			std::vector<SphereData_Spectral*> &io_prog_fields,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	);


	void run_timestep_1(
			SphereData_Spectral &io_prognostic_field,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	);


	void run_timestep_2(
			std::vector<SphereData_Spectral*> &io_prog_fields,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	);



	void run_timestep_3(
			std::vector<SphereData_Spectral*> &io_prog_fields,	///< prognostic variables
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_fixed_dt,
			double i_simulation_timestamp,

			// for varying velocity fields
			const BenchmarksSphereAdvection *i_sphereBenchmarks
	);



	virtual ~SphereAdvection_TS_na_sl();
};

#endif
