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

		double i_dt,						///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const SWESphereBenchmarksCombined *i_sphereBenchmarks,
		SphereData_Physical &io_U_phi_phys
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	SphereData_Spectral &U_phi = io_U_phi;
	SphereData_Spectral &U_vrt = io_U_vrt;
	SphereData_Spectral &U_div = io_U_div;

	if (i_simulation_timestamp == 0)
	{
		U_phi_prev = U_phi;
		U_vrt_prev = U_vrt;
		U_div_prev = U_div;
	}

#if 1
	// shouldn't be necessary
	if (i_sphereBenchmarks != nullptr)
	{
		i_sphereBenchmarks->update_time_varying_fields_pert(U_phi, U_vrt, U_div, i_simulation_timestamp);
		i_sphereBenchmarks->update_time_varying_fields_pert(U_phi_prev, U_vrt_prev, U_div_prev, i_simulation_timestamp - i_dt);
	}
#endif

	// IMPORTANT!!! WE DO NOT USE THE ROBERT TRANSFORMATION HERE!!!

	SphereData_Physical U_u(sphereDataConfig);
	SphereData_Physical U_v(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt, U_div, U_u, U_v, false);

	SphereData_Physical U_u_prev(sphereDataConfig);
	SphereData_Physical U_v_prev(sphereDataConfig);
	op.vortdiv_to_uv(U_vrt_prev, U_div_prev, U_u_prev, U_v_prev, false);

	// OUTPUT: position of departure points at t
	ScalarDataArray pos_lon_d(sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray pos_lat_d(sphereDataConfig->physical_array_data_number_of_elements);

#if 0

	int num_elements = sphereDataConfig->physical_array_data_number_of_elements;

	/*
	 * Compute Cartesian velocity
	 */
	ScalarDataArray u_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(U_u);
	ScalarDataArray v_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(U_v);

	ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
	SWEETMath::latlon_velocity_to_cartesian_velocity(
		semiLagrangian.pos_lon_A, semiLagrangian.pos_lat_A,
		u_lon_array, v_lat_array,
		&vel_x_A, &vel_y_A, &vel_z_A
	);

	/*
	 * Do advection in Cartesian space
	 */
	double dt_div_radius = i_dt / simVars.sim.sphere_radius;

	ScalarDataArray new_pos_x_d(num_elements), new_pos_y_d(num_elements), new_pos_z_d(num_elements);
	semiLagrangian.doAdvectionOnSphere(
		semiLagrangian.pos_x_A, semiLagrangian.pos_y_A, semiLagrangian.pos_z_A,
		-dt_div_radius*vel_x_A, -dt_div_radius*vel_y_A, -dt_div_radius*vel_z_A,
		new_pos_x_d, new_pos_y_d, new_pos_z_d,

		false
	);

	/*
	 * Departure point to lat/lon coordinate
	 */
	SWEETMath::cartesian_to_latlon(
			new_pos_x_d, new_pos_y_d, new_pos_z_d,
			pos_lon_d, pos_lat_d
		);

	SphereData_Physical new_prog_phi_phys(sphereDataConfig);
	sampler2D.bicubic_scalar(
			io_U_phi.getSphereDataPhysical(),
			pos_lon_d, pos_lat_d,
			new_prog_phi_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			false
	);

	io_U_phi = new_prog_phi_phys;
	return;

#endif

	//double dt_div_radius = simVars.timecontrol.current_timestep_size / simVars.sim.sphere_radius;
	double one_div_radius = 1.0 / simVars.sim.sphere_radius;

	semiLagrangian.semi_lag_departure_points_settls(
			//dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			//dt_div_radius*U_u, dt_div_radius*U_v,
			U_u_prev, U_v_prev,
			U_u, U_v,

			simVars.timecontrol.current_timestep_size,
			simVars.timecontrol.current_simulation_time,
			simVars.sim.sphere_radius,
			i_sphereBenchmarks,

			pos_lon_d, pos_lat_d,

			op,

			timestepping_order,
			simVars.disc.semi_lagrangian_max_iterations,
			simVars.disc.semi_lagrangian_convergence_threshold,
			simVars.disc.semi_lagrangian_approximate_sphere_geometry
	);

	U_phi_prev = U_phi;
	U_vrt_prev = U_vrt;
	U_div_prev = U_div;

	SphereData_Physical new_prog_phi_physx(sphereDataConfig);

	sampler2D.bicubic_scalar(
			io_U_phi.getSphereDataPhysical(),
			pos_lon_d,
			pos_lat_d,
			new_prog_phi_physx,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			false
	);

	io_U_phi = new_prog_phi_physx;
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

