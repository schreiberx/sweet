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




void SphereAdvection_TS_na_sl::interpolate_departure_point_uvw(
		const SphereData_Spectral &i_u,
		const SphereData_Spectral &i_v,
		const SphereData_Spectral &i_w,

		const ScalarDataArray &i_pos_lon_D,
		const ScalarDataArray &i_pos_lat_D,

		SphereData_Spectral &o_u,
		SphereData_Spectral &o_v,
		SphereData_Spectral &o_w
)
{
	const SphereData_Config *sphereDataConfig = i_u.sphereDataConfig;

	o_u.setup_if_required(i_u.sphereDataConfig);
	o_v.setup_if_required(i_v.sphereDataConfig);
	o_w.setup_if_required(i_w.sphereDataConfig);

	/*
	 * First we sample the field at the departure point
	 */
	SphereData_Physical u_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			o_u.toPhys(),
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical v_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			o_v.toPhys(),
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	SphereData_Physical w_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
			o_w.toPhys(),
			i_pos_lon_D, i_pos_lat_D,
			true,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
		);

	/*
	 * Convert to scalar data arrays
	 */
	ScalarDataArray V_u_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(u_tmp_D);
	ScalarDataArray V_v_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(v_tmp_D);
	ScalarDataArray V_w_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(w_tmp_D);



	/*
	 * Here we have the velocity at the departure points.
	 *
	 * Now we need to take the rotation of this vector into account!
	 */

	/*
	 * Compute departure position
	 */
	ScalarDataArray P_x_D, P_y_D, P_z_D;
	SWEETMath::latlon_to_cartesian(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

	ScalarDataArray	&P_x_A = semiLagrangian.pos_x_A;
	ScalarDataArray	&P_y_A = semiLagrangian.pos_y_A;
	ScalarDataArray	&P_z_A = semiLagrangian.pos_z_A;

	/*
	 * Compute rotation angle based on departure and arrival position
	 */
	ScalarDataArray rotation_angle_ =
			SWEETMath::dot_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A
			);

	// Can be slightly larger than 1 due to round-off issues, leading to NaN, hence this hack
	rotation_angle_ = SWEETMath::min(rotation_angle_, 1.0);

	/*
	 * Compute rotation angle
	 */
	ScalarDataArray rotation_angle = SWEETMath::arccos(rotation_angle_);

	/*
	 * Compute Rotation axis and normalize
	 */
	ScalarDataArray rot_x, rot_y, rot_z;
	SWEETMath::cross_prod(
			P_x_D, P_y_D, P_z_D,
			P_x_A, P_y_A, P_z_A,
			rot_x, rot_y, rot_z
		);
	SWEETMath::normalize_with_threshold(rot_x, rot_y, rot_z);

	/*
	 * Rotate vector (using transpose of rotation matrix without translation!)
	 */
	ScalarDataArray V_u_A, V_v_A, V_w_A;
	SWEETMath::rotate_3d_vector_normalized_rotation_axis(
			V_u_D, V_v_D, V_w_D,
			rotation_angle,
			rot_x, rot_y, rot_z,
			V_u_A, V_v_A, V_w_A
		);

	o_u = Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_u_A, sphereDataConfig);
	o_v = Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_v_A, sphereDataConfig);
	o_w = Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_w_A, sphereDataConfig);
}


void SphereAdvection_TS_na_sl::run_timestep(
		std::vector<SphereData_Spectral*> &io_prognostic_fields,	///< prognostic variables
		SphereData_Physical &io_U_u,
		SphereData_Physical &io_U_v,

		double i_dt,
		double i_simulation_timestamp,

		// for varying velocity fields
		const BenchmarksSphereAdvection *i_sphereBenchmarks
)
{
	if (io_prognostic_fields.size() == 1)
	{
		run_timestep(
				*io_prognostic_fields[0],
				io_U_u,
				io_U_v,
				i_dt,
				i_simulation_timestamp,
				i_sphereBenchmarks
			);
		return;
	}

	SWEETAssert(io_prognostic_fields.size() == 3, "Only 3D Vector supported");

	const SphereData_Config *sphereDataConfig = io_prognostic_fields[0]->sphereDataConfig;

	if (i_simulation_timestamp == 0)
	{
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

	// OUTPUT: position of departure points at t
	ScalarDataArray pos_lon_d(sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray pos_lat_d(sphereDataConfig->physical_array_data_number_of_elements);

	double dt_div_radius = simVars.timecontrol.current_timestep_size / simVars.sim.sphere_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_d, pos_lat_d
	);



	interpolate_departure_point_uvw(
			*io_prognostic_fields[0],
			*io_prognostic_fields[1],
			*io_prognostic_fields[2],

			pos_lon_d,
			pos_lat_d,

			*io_prognostic_fields[0],
			*io_prognostic_fields[1],
			*io_prognostic_fields[2]
	);
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

	if (i_simulation_timestamp == 0)
	{
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

	// OUTPUT: position of departure points at t
	ScalarDataArray pos_lon_D(sphereDataConfig->physical_array_data_number_of_elements);
	ScalarDataArray pos_lat_D(sphereDataConfig->physical_array_data_number_of_elements);

	double dt_div_radius = simVars.timecontrol.current_timestep_size / simVars.sim.sphere_radius;

	semiLagrangian.semi_lag_departure_points_settls_specialized(
			dt_div_radius*U_u_prev, dt_div_radius*U_v_prev,
			dt_div_radius*io_U_u, dt_div_radius*io_U_v,

			pos_lon_D, pos_lat_D
	);

	U_u_prev = io_U_u;
	U_v_prev = io_U_v;

	SphereData_Physical new_prog_phi_phys =
		sphereSampler.bicubic_scalar_ret_phys(
			io_U_phi.toPhys(),
			pos_lon_D,
			pos_lat_D,
			false,
			simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points,
			simVars.disc.semi_lagrangian_interpolation_limiter
	);

	io_U_phi = new_prog_phi_phys;
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

