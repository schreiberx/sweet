/*
 * Adv_Sphere_TS_na_sl.cpp
 *
 *  Created on: 29 Mar 2018
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "Adv_Sphere_TS_na_sl.hpp"




void Adv_Sphere_TS_na_sl::run_timestep(
		SphereData_Spectral &io_U_phi,		///< prognostic variables
		SphereData_Spectral &io_U_vrt,		///< prognostic variables
		SphereData_Spectral &io_U_div,		///< prognostic variables

		double i_fixed_dt,		///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const SWESphereBenchmarksCombined *i_sphereBenchmarks
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	SphereData_Spectral &U_phi = io_U_phi;
	SphereData_Spectral &U_vrt = io_U_vrt;
	SphereData_Spectral &U_div = io_U_div;

	if (i_fixed_dt <= 0)
		FatalError("Only constant time step size allowed");

	double dt = simVars.timecontrol.current_timestep_size;

	if (i_simulation_timestamp == 0)
	{
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}
	else
	{
#if 0
		// shouldn't be necessary
		if (i_sphereBenchmarks != nullptr)
		{
			i_sphereBenchmarks->update_time_varying_fields_pert(U_phi_prev, U_vrt_prev, U_div_prev, i_simulation_timestamp - i_fixed_dt);
		}
#endif
	}

	// IMPORTANT!!! WE DO NOT USE THE ROBERT TRANSFORMATION HERE!!!

	SphereData_Physical U_u(sphereDataConfig);
	SphereData_Physical U_v(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt, U_div, U_u, U_v, false);

	SphereData_Physical U_u_prev(sphereDataConfig);
	SphereData_Physical U_v_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt_prev, U_div_prev, U_u_prev, U_v_prev, false);

	// OUTPUT: position of departure points at t
	ScalarDataArray posx_d(io_U_phi.sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray posy_d(io_U_phi.sphereDataConfig->physical_array_data_number_of_elements);

	double dt_div_radius = simVars.timecontrol.current_timestep_size / simVars.sim.sphere_radius;

	semiLagrangian.semi_lag_departure_points_settls(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*U_u, dt_div_radius*U_v,

			posx_a, posy_a,
			posx_d, posy_d,

			timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry
	);

	U_phi_prev = U_phi;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	SphereData_Physical new_prog_phi_phys(sphereDataConfig);

	sampler2D.bicubic_scalar(
			io_U_phi.getSphereDataPhysical(),
			posx_d,
			posy_d,
			new_prog_phi_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
	);

	io_U_phi = new_prog_phi_phys;
}



/*
 * Setup
 */
void Adv_Sphere_TS_na_sl::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		FatalError("Only 1st and 2nd order for SL integration supported");

	const SphereData_Config *sphereDataConfig = op.sphereDataConfig;

	posx_a.setup(sphereDataConfig->physical_array_data_number_of_elements);
	posy_a.setup(sphereDataConfig->physical_array_data_number_of_elements);

	// setup some test sampling points
	// we use 2 arrays - one for each sampling position

	posx_a.update_lambda_array_indices(
		[&](int idx, double &io_data)
		{
			int i = idx % sphereDataConfig->physical_num_lon;
			//int j = idx / sphereDataConfig->physical_data_size[0];

			io_data = 2.0*M_PI*(double)i/(double)sphereDataConfig->physical_num_lon;
			assert(io_data >= 0);
			assert(io_data < 2.0*M_PI);
		}
	);
	posy_a.update_lambda_array_indices(
			[&](int idx, double &io_data)
		{
			//int i = idx % sphereDataConfig->physical_data_size[0];
			int j = idx / sphereDataConfig->physical_num_lon;

			io_data = sphereDataConfig->lat[j];

			assert(io_data >= -M_PI*0.5);
			assert(io_data <= M_PI*0.5);
		}
	);

	sampler2D.setup(sphereDataConfig);

	//PXT- This just calls sampler2D.setup, so any reason for having it?
	semiLagrangian.setup(sphereDataConfig, simVars);
}


Adv_Sphere_TS_na_sl::Adv_Sphere_TS_na_sl(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op)
{
	setup(simVars.disc.timestepping_order);
}



Adv_Sphere_TS_na_sl::~Adv_Sphere_TS_na_sl()
{
}

