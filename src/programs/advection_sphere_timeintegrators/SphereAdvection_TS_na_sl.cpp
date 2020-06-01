/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 */

#include "../advection_sphere_timeintegrators/SphereAdvection_TS_na_sl.hpp"



bool SphereAdvection_TS_na_sl::implements_timestepping_method(const std::string &i_timestepping_method)
{
	return i_timestepping_method == "na_sl";
}

std::string SphereAdvection_TS_na_sl::string_id()
{
	return "na_sl";
}



void SphereAdvection_TS_na_sl::setup_auto()
{
	setup(simVars.disc.timestepping_order);
}

std::string SphereAdvection_TS_na_sl::get_help()
{
	std::ostringstream stream;
	stream << " + SphereAdvection_TS_na_sl:" << std::endl;
	stream << "    * 'na_sl'" << std::endl;

	return stream.str();
}


void SphereAdvection_TS_na_sl::run_timestep(
		std::vector<SphereData_Spectral*> &io_prog_fields,		///< prognostic variables
		SphereData_Physical &io_U_u,
		SphereData_Physical &io_U_v,

		double i_dt,						///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const BenchmarksSphereAdvection *i_sphereBenchmarks
)
{

}

void SphereAdvection_TS_na_sl::run_timestep(
		SphereData_Spectral &io_U_phi,		///< prognostic variables
		SphereData_Physical &io_U_u,
		SphereData_Physical &io_U_v,

		double i_dt,						///< if this value is not equal to 0, use this time step size instead of computing one
		double i_simulation_timestamp,

		// for varying velocity fields
		const BenchmarksSphereAdvection *i_sphereBenchmarks
)
{
	const SphereData_Config *sphereDataConfig = io_U_phi.sphereDataConfig;

	SphereData_Spectral &U_phi = io_U_phi;

	if (i_simulation_timestamp == 0)
	{
		U_phi_prev = U_phi;
		U_u_prev = io_U_u;
		U_v_prev = io_U_v;
	}

	/*
	 * For time-varying fields, update the vrt/div field based on the given simulation timestamp
	 */
	if (i_sphereBenchmarks)
	{
		i_sphereBenchmarks->master->get_varying_velocities(io_U_u, io_U_v, i_simulation_timestamp);
		i_sphereBenchmarks->master->get_varying_velocities(U_u_prev, U_v_prev, i_simulation_timestamp - i_dt);
	}


	// IMPORTANT!!! WE DO NOT USE THE ROBERT TRANSFORMATION HERE!!!

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
		vel_x_A, vel_y_A, vel_z_A
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
	sphereSampler.bicubic_scalar(
			io_U_phi.toPhys(),
			pos_lon_d, pos_lat_d,
			new_prog_phi_phys,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter,
			false
	);

	io_U_phi = new_prog_phi_phys;
	return;

#endif

	double dt_div_radius = simVars.timecontrol.current_timestep_size / simVars.sim.sphere_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_d, pos_lat_d
	);

	U_phi_prev = U_phi;
	U_u_prev = io_U_u;
	U_v_prev = io_U_v;

	SphereData_Physical new_prog_phi_physx =
		sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi.toPhys(),
			pos_lon_d,
			pos_lat_d,
			false,
			false,
			simVars.disc.semi_lagrangian_interpolation_limiter
	);

	io_U_phi = new_prog_phi_physx;
}



/*
 * Setup
 */
void SphereAdvection_TS_na_sl::setup(
		int i_order	///< order of RK time stepping method
)
{
	timestepping_order = i_order;

	if (timestepping_order > 2 || timestepping_order <= 0)
		SWEETError("Only 1st and 2nd order for SL integration supported");
}


SphereAdvection_TS_na_sl::SphereAdvection_TS_na_sl(
		SimulationVariables &i_simVars,
		SphereOperators_SphereData &i_op
)	:
		simVars(i_simVars),
		op(i_op),
		semiLagrangian(simVars),
		sphereSampler(semiLagrangian.sphereSampler)
{
	setup(simVars.disc.timestepping_order);

	semiLagrangian.setup(op.sphereDataConfig);
}



SphereAdvection_TS_na_sl::~SphereAdvection_TS_na_sl()
{
}

