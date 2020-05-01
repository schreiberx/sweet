/*
 * SphereDataSemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 *  Updated to sphere on 28th March 2018
 */
#ifndef SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_
#define SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_

#include <limits>
#include <sweet/ScalarDataArray.hpp>
#include <sweet/SWEETMath.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/Convert_ScalarDataArray_to_SphereDataPhysical.hpp>
#include <sweet/SimulationVariables.hpp>

#include <sweet/sphere/Convert_SphereDataPhysical_to_ScalarDataArray.hpp>
#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>



class SphereTimestepping_SemiLagrangian
{
	SimulationVariables &simVars;

	const SphereData_Config *sphereDataConfig;

	SphereData_Physical sl_coriolis;

	// Arrival points

public:
	ScalarDataArray pos_lat_A;
	ScalarDataArray pos_lon_A;

	ScalarDataArray pos_x_A;
	ScalarDataArray pos_y_A;
	ScalarDataArray pos_z_A;

	SphereOperators_Sampler_SphereDataPhysical sphereSampler;


	int timestepping_order;
	int semi_lagrangian_max_iterations;
	double semi_lagrangian_convergence_threshold;
	bool semi_lagrangian_approximate_sphere_geometry;
	bool semi_lagrangian_interpolation_limiter;

	enum EnumTrajectories
	{
		E_TRAJECTORY_METHOD_CANONICAL,
		E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE,
		E_TRAJECTORY_METHOD_SETTLS_HORTAL
	};


	EnumTrajectories trajectory_method;


	SphereTimestepping_SemiLagrangian(
			SimulationVariables &i_simVars
	)	:
		simVars(i_simVars),
		sphereDataConfig(nullptr)
	{
	}


	SphereTimestepping_SemiLagrangian(
			SimulationVariables &i_simVars,
			const SphereData_Config *i_sphereDataConfig
	)	:
		simVars(i_simVars)
	{
		setup(i_sphereDataConfig);
	}


	void setup(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereSampler.setup(sphereDataConfig);

		if (simVars.disc.semi_lagrangian_departure_point_method == "canonical")
		{
			trajectory_method = E_TRAJECTORY_METHOD_CANONICAL;
		}
		else if (simVars.disc.semi_lagrangian_departure_point_method == "midpoint")
		{
			trajectory_method = E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE;
		}
		else if (simVars.disc.semi_lagrangian_departure_point_method == "settls")
		{
			trajectory_method = E_TRAJECTORY_METHOD_SETTLS_HORTAL;
		}
		else
		{
			FatalError(std::string("Trajectory method '")+simVars.disc.semi_lagrangian_departure_point_method+"' not supported");
		}


		timestepping_order = simVars.disc.timestepping_order;
		semi_lagrangian_max_iterations = simVars.disc.semi_lagrangian_max_iterations;
		semi_lagrangian_convergence_threshold = simVars.disc.semi_lagrangian_convergence_threshold;
		semi_lagrangian_approximate_sphere_geometry = simVars.disc.semi_lagrangian_approximate_sphere_geometry;
		semi_lagrangian_interpolation_limiter = simVars.disc.semi_lagrangian_interpolation_limiter;



		pos_lon_A.setup(sphereDataConfig->physical_array_data_number_of_elements);
		pos_lon_A.update_lambda_array_indices(
			[&](int idx, double &io_data)
			{
				int i = idx % sphereDataConfig->physical_num_lon;

				io_data = 2.0*M_PI*(double)i/(double)sphereDataConfig->physical_num_lon;
				assert(io_data >= 0);
				assert(io_data < 2.0*M_PI);
			}
		);

		pos_lat_A.setup(sphereDataConfig->physical_array_data_number_of_elements);
		pos_lat_A.update_lambda_array_indices(
				[&](int idx, double &io_data)
			{
				int j = idx / sphereDataConfig->physical_num_lon;

				io_data = sphereDataConfig->lat[j];

				assert(io_data >= -M_PI*0.5);
				assert(io_data <= M_PI*0.5);
			}
		);

		pos_x_A.setup(i_sphereDataConfig->physical_array_data_number_of_elements);
		pos_y_A.setup(i_sphereDataConfig->physical_array_data_number_of_elements);
		pos_z_A.setup(i_sphereDataConfig->physical_array_data_number_of_elements);

		SWEETMath::latlon_to_cartesian(
				pos_lon_A, pos_lat_A,
				pos_x_A, pos_y_A, pos_z_A
			);


		//SphereData_Physical u_lon_
		if (simVars.sim.sphere_use_fsphere)
		{
			FatalError("Not supported");
		}


		sl_coriolis.setup_if_required(sphereDataConfig);
		sl_coriolis.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*simVars.sim.sphere_rotating_coriolis_omega*simVars.sim.sphere_radius;
			}
		);
	}



	/**
	 * Do 1st order accurate advection on the sphere for given
	 *
	 *  - starting point,
	 *  - velocity*dt
	 *
	 * in Cartesian space 
	 */
	inline static
	void doAdvectionOnSphere(
		const ScalarDataArray &i_pos_x,
		const ScalarDataArray &i_pos_y,
		const ScalarDataArray &i_pos_z,

		const ScalarDataArray &i_dt_velocity_x,
		const ScalarDataArray &i_dt_velocity_y,
		const ScalarDataArray &i_dt_velocity_z,

		ScalarDataArray &o_pos_x,
		ScalarDataArray &o_pos_y,
		ScalarDataArray &o_pos_z,

		double i_approximate_sphere_geometry
	)
	{
		if (i_approximate_sphere_geometry)
		{
			/*
			 * This just uses an approximation of the sphere geometry.
			 */

			/*
			 * Step 1) Apply the velocity vector in Cartesian space, ignoring the sphere's curvature
			 */
			o_pos_x = i_pos_x + i_dt_velocity_x;
			o_pos_y = i_pos_y + i_dt_velocity_y;
			o_pos_z = i_pos_z + i_dt_velocity_z;

			/*
			 * Step 2) Normalize the resulting position
			 */
			SWEETMath::normalize(o_pos_x, o_pos_y, o_pos_z);
			return;
		}


		/*
		 * This version implements an accurate geometry.
		 */
		{
			/*
			 * Step 1) Compute rotation axis with cross product
			 */
			ScalarDataArray rotation_axis_x(i_pos_x.number_of_elements);
			ScalarDataArray rotation_axis_y(i_pos_x.number_of_elements);
			ScalarDataArray rotation_axis_z(i_pos_x.number_of_elements);

			SWEETMath::cross_prod(
					i_pos_x, i_pos_y, i_pos_z,
					i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z,
					rotation_axis_x, rotation_axis_y, rotation_axis_z
				);

			/*
			 * Normalize rotation axis since it's likely not normalized yet
			 */
			SWEETMath::normalize_threshold(
					rotation_axis_x,
					rotation_axis_y,
					rotation_axis_z,
					1e-12				// close-to-0 threshold
				);

			/*
			 * Compute rotation angle
			 *
			 * No rescaling by 1/(2pi) since the change of angle is directly
			 * given by the magnitude of the angular velocity by its definition
			 */
			ScalarDataArray angle = SWEETMath::length(i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z);

			/*
			 * Rotate
			 */
			SWEETMath::rotate_3d_normalized_rotation_axis(
					i_pos_x,
					i_pos_y,
					i_pos_z,
					angle,
					rotation_axis_x,
					rotation_axis_y,
					rotation_axis_z,
					o_pos_x,
					o_pos_y,
					o_pos_z
				);
		}
	}

#if 0
	/*
	 * Compute SL departure points on unit sphere for given dt*(u,v) velocities
	 */
	void semi_lag_departure_points_settls(
		const SphereData_Physical &i_u_lon_prev,	///< Velocities at time t-1
		const SphereData_Physical &i_v_lat_prev,

		const SphereData_Physical &i_dt_u_lon, 		///< Velocities at time t
		const SphereData_Physical &i_dt_v_lat,

		double i_dt,
		double i_timestamp,
		double i_radius,

		// for varying velocity fields
		const SWESphereBenchmarksCombined *i_sphereBenchmarks,

		ScalarDataArray &o_pos_lon_D, 	///< OUTPUT: Position of departure points x / y
		ScalarDataArray &o_pos_lat_D,

		SphereOperators_SphereData &op,

		int i_timestepping_order,

		int max_iterations = 2,
		double i_convergence_tolerance = -1,
		int i_approximate_sphere_geometry = 0,
		bool use_interpolation_limiters = false
	)
	{
		o_pos_lon_D.setup_if_required(pos_lon_A);
		o_pos_lat_D.setup_if_required(pos_lon_A);

		std::size_t num_elements = o_pos_lon_D.number_of_elements;

		if (i_timestepping_order == 1)
		{
			/*
			 * Compute Cartesian velocity
			 */
			ScalarDataArray u_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
			ScalarDataArray v_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

			ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
			SWEETMath::latlon_velocity_to_cartesian_velocity(
					pos_lon_A, pos_lat_A,
					u_lon_array, v_lat_array,
					vel_x_A, vel_y_A, vel_z_A
				);

			/*
			 * Do advection in Cartesian space
			 */

			double dt_div_radius = i_dt/i_radius;
			ScalarDataArray new_pos_x_d(num_elements), new_pos_y_d(num_elements), new_pos_z_d(num_elements);
			doAdvectionOnSphere(
				pos_x_A, pos_y_A, pos_z_A,
				-dt_div_radius*vel_x_A, -dt_div_radius*vel_y_A, -dt_div_radius*vel_z_A,
				new_pos_x_d, new_pos_y_d, new_pos_z_d,

				i_approximate_sphere_geometry
			);

			/*
			 * Departure point to lat/lon coordinate
			 */
			SWEETMath::cartesian_to_latlon(
					new_pos_x_d, new_pos_y_d, new_pos_z_d,
					o_pos_lon_D, o_pos_lat_D
				);
			return;
		}


		if (i_timestepping_order == 2)
		{
			if (trajectory_method == E_TRAJECTORY_METHOD_CANONICAL)
			{
				/*
				 * Standard iterative method
				 *
				 * See also Michail Diamantarkis paper, p. 185
				 */

				/*
				 * Prepare
				 */
				const SphereData_Physical &vel_lon = i_dt_u_lon;
				const SphereData_Physical &vel_lat = i_dt_v_lat;

				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(vel_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(vel_lat);

				/*
				 * Polar => Cartesian velocities
				 */
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Step 1)
				 */
				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				double dt_div_radius = i_dt/i_radius;
				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-dt_div_radius*vel_x_A,
					-dt_div_radius*vel_y_A,
					-dt_div_radius*vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					i_approximate_sphere_geometry
				);


				/*
				 * 2 iterations to get midpoint
				 */
				for (int i = 0; i < 2; i++)
				{
					/*
					 * Step 2a
					 */
					ScalarDataArray pos_x_mid = 0.5*(pos_x_A + pos_x_D);
					ScalarDataArray pos_y_mid = 0.5*(pos_y_A + pos_y_D);
					ScalarDataArray pos_z_mid = 0.5*(pos_z_A + pos_z_D);

					ScalarDataArray pos_lon_mid(sphereDataConfig->physical_array_data_number_of_elements);
					ScalarDataArray pos_lat_mid(sphereDataConfig->physical_array_data_number_of_elements);

					SWEETMath::cartesian_to_latlon(
							pos_x_mid, pos_y_mid, pos_z_mid,
							pos_lon_mid, pos_lat_mid
						);

					ScalarDataArray vel_u_mid = sphereSampler.bilinear_scalar(vel_lon, pos_lon_mid, pos_lat_mid, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray vel_v_mid = sphereSampler.bilinear_scalar(vel_lat, pos_lon_mid, pos_lat_mid, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_mid(num_elements), vel_y_mid(num_elements), vel_z_mid(num_elements);
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							pos_lon_mid, pos_lat_mid,
							vel_u_mid, vel_v_mid,
							vel_x_mid, vel_y_mid, vel_z_mid
						);

					/*
					 * Step 2b
					 */
					double dt_div_radius = i_dt/i_radius;
					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-dt_div_radius*vel_x_mid,
						-dt_div_radius*vel_y_mid,
						-dt_div_radius*vel_z_mid,

						pos_x_D,
						pos_y_D,
						pos_z_D,

						i_approximate_sphere_geometry
					);
				}

				// convert final points from Cartesian space to angular space
				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE)
			{

				/*
				 * Ritchies midpoint rule
				 */
				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */
				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				double dt_div_radius = i_dt/i_radius;
				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-0.5*dt_div_radius*vel_x_A,
					-0.5*dt_div_radius*vel_y_A,
					-0.5*dt_div_radius*vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					i_approximate_sphere_geometry
				);


				/*
				 * 2 iterations to get midpoint
				 */
				for (int i = 0; i < 2; i++)
				{
					SWEETMath::cartesian_to_latlon(
							pos_x_D, pos_y_D, pos_z_D,
							o_pos_lon_D, o_pos_lat_D
						);

					ScalarDataArray u_D = sphereSampler.bilinear_scalar(i_dt_u_lon, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray v_D = sphereSampler.bilinear_scalar(i_dt_v_lat, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_D(num_elements), vel_y_D(num_elements), vel_z_D(num_elements);
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							o_pos_lon_D, o_pos_lat_D,
							u_D, v_D,
							vel_x_D, vel_y_D, vel_z_D
						);

					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-0.5*i_dt*vel_x_D,
						-0.5*i_dt*vel_y_D,
						-0.5*i_dt*vel_z_D,

						pos_x_D,
						pos_y_D,
						pos_z_D,

						i_approximate_sphere_geometry
					);
				}

				/*
				 * Given the midpoint at pos_?_D, we compute the full time step
				 */

				ScalarDataArray dot2 = 2.0*(pos_x_D*pos_x_A + pos_y_D*pos_y_A + pos_z_D*pos_z_A);

				pos_x_D = dot2*pos_x_D - pos_x_A;
				pos_y_D = dot2*pos_y_D - pos_y_A;
				pos_z_D = dot2*pos_z_D - pos_z_A;

				// convert final points from Cartesian space to angular space
				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_SETTLS_HORTAL)
			{
				/**
				 * See SETTLS paper
				 * Hortal, M. (2002). The development and testing of a new two-time-level semi-Lagrangian scheme (SETTLS) in the ECMWF forecast model. Q. J. R. Meteorol. Soc., 2, 1671–1687.
				 *
				 * We use the SETTLS formulation also to compute the departure points.
				 */

				// Extrapolate velocities at departure points
				SphereData_Physical u_extrapol = 2.0*i_dt_u_lon - i_u_lon_prev;
				SphereData_Physical v_extrapol = 2.0*i_dt_v_lat - i_v_lat_prev;

				// Convert velocities along lon/lat to scalardata array
				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */

				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				double dt_div_radius = i_dt/i_radius;
				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-dt_div_radius*vel_x_A,
					-dt_div_radius*vel_y_A,
					-dt_div_radius*vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					i_approximate_sphere_geometry
				);

				double diff = 999999;

				for (int iters = 0; iters < max_iterations; iters++)
				{
					SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);


					/*
					 * WARNING: Never convert this to vort/div space!!!
					 * This creates some artificial waves
					 */
					ScalarDataArray vel_lon_extrapol_D = sphereSampler.bilinear_scalar(u_extrapol, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray vel_lat_extrapol_D = sphereSampler.bilinear_scalar(v_extrapol, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_extrapol_D(num_elements), vel_y_extrapol_D(num_elements), vel_z_extrapol_D(num_elements);

					// polar => Cartesian coordinates
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							o_pos_lon_D, o_pos_lat_D,
							vel_lon_extrapol_D, vel_lat_extrapol_D,
							vel_x_extrapol_D, vel_y_extrapol_D, vel_z_extrapol_D
						);

					ScalarDataArray new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					double dt_div_radius = i_dt/i_radius;
					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-dt_div_radius*0.5*(vel_x_extrapol_D + vel_x_A),
						-dt_div_radius*0.5*(vel_y_extrapol_D + vel_y_A),
						-dt_div_radius*0.5*(vel_z_extrapol_D + vel_z_A),

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						i_approximate_sphere_geometry
					);

					if (i_convergence_tolerance > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
								(pos_y_D-new_pos_y_D).reduce_maxAbs() +
								(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < i_convergence_tolerance)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (i_convergence_tolerance > 0)
				{
					if (diff > i_convergence_tolerance)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << i_convergence_tolerance << std::endl;
					}
				}

				// convert final points from Cartesian space to angular space
				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else
			{
				FatalError("Unknown departure point calculation method");
			}

			return;
		}

		FatalError("Only 1st and 2nd order time integration supported");
	}
#endif


	/*
	 * Compute SL departure points on unit sphere for given dt*(u,v) velocities
	 */
	void semi_lag_departure_points_settls_specialized(
		const SphereData_Physical &i_u_lon_prev,	///< Velocities at time t-1
		const SphereData_Physical &i_v_lat_prev,

		const SphereData_Physical &i_dt_u_lon, 		///< Velocities at time t
		const SphereData_Physical &i_dt_v_lat,

		ScalarDataArray &o_pos_lon_D, 	///< OUTPUT: Position of departure points x / y
		ScalarDataArray &o_pos_lat_D
	)
	{
		o_pos_lon_D.setup_if_required(pos_lon_A);
		o_pos_lat_D.setup_if_required(pos_lon_A);

		std::size_t num_elements = o_pos_lon_D.number_of_elements;

		if (timestepping_order == 1)
		{
			/*
			 * Compute Cartesian velocity
			 */
			ScalarDataArray u_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
			ScalarDataArray v_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

			ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
			SWEETMath::latlon_velocity_to_cartesian_velocity(
					pos_lon_A, pos_lat_A,
					u_lon_array, v_lat_array,
					vel_x_A, vel_y_A, vel_z_A
				);

			/*
			 * Do advection in Cartesian space
			 */

			ScalarDataArray new_pos_x_d(num_elements), new_pos_y_d(num_elements), new_pos_z_d(num_elements);
			doAdvectionOnSphere(
				pos_x_A, pos_y_A, pos_z_A,
				-vel_x_A, -vel_y_A, -vel_z_A,
				new_pos_x_d, new_pos_y_d, new_pos_z_d,

				semi_lagrangian_approximate_sphere_geometry
			);

			/*
			 * Departure point to lat/lon coordinate
			 */
			SWEETMath::cartesian_to_latlon(
					new_pos_x_d, new_pos_y_d, new_pos_z_d,
					o_pos_lon_D, o_pos_lat_D
				);
			return;
		}


		if (timestepping_order == 2)
		{
			if (trajectory_method == E_TRAJECTORY_METHOD_CANONICAL)
			{
				/*
				 * Standard iterative method
				 *
				 * See also Michail Diamantarkis paper, p. 185
				 */

				/*
				 * Prepare
				 */
				const SphereData_Physical &vel_lon = i_dt_u_lon;
				const SphereData_Physical &vel_lat = i_dt_v_lat;

				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(vel_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(vel_lat);

				/*
				 * Polar => Cartesian velocities
				 */
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Step 1)
				 */
				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-vel_x_A,
					-vel_y_A,
					-vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere_geometry
				);


				/*
				 * 2 iterations to get midpoint
				 */
				double diff = -1;
				for (int i = 0; i < semi_lagrangian_max_iterations; i++)
				{
					/*
					 * Step 2a
					 */
					ScalarDataArray pos_x_mid = 0.5*(pos_x_A + pos_x_D);
					ScalarDataArray pos_y_mid = 0.5*(pos_y_A + pos_y_D);
					ScalarDataArray pos_z_mid = 0.5*(pos_z_A + pos_z_D);

					ScalarDataArray pos_lon_mid(sphereDataConfig->physical_array_data_number_of_elements);
					ScalarDataArray pos_lat_mid(sphereDataConfig->physical_array_data_number_of_elements);

					SWEETMath::cartesian_to_latlon(
							pos_x_mid, pos_y_mid, pos_z_mid,
							pos_lon_mid, pos_lat_mid
						);

					ScalarDataArray vel_u_mid = sphereSampler.bilinear_scalar(vel_lon, pos_lon_mid, pos_lat_mid, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray vel_v_mid = sphereSampler.bilinear_scalar(vel_lat, pos_lon_mid, pos_lat_mid, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_mid(num_elements), vel_y_mid(num_elements), vel_z_mid(num_elements);
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							pos_lon_mid, pos_lat_mid,
							vel_u_mid, vel_v_mid,
							vel_x_mid, vel_y_mid, vel_z_mid
						);

					// convert final points from Cartesian space to angular space
					ScalarDataArray new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					/*
					 * Step 2b
					 */
					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-vel_x_mid,
						-vel_y_mid,
						-vel_z_mid,

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere_geometry
					);


					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
										(pos_y_D-new_pos_y_D).reduce_maxAbs() +
										(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				// convert final points from Cartesian space to angular space
				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE)
			{

				/*
				 * Ritchies midpoint rule
				 */
				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */
				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-0.5*vel_x_A,
					-0.5*vel_y_A,
					-0.5*vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere_geometry
				);


				/*
				 * 2 iterations to get midpoint
				 */
				double diff = -1;
				for (int i = 0; i < semi_lagrangian_max_iterations; i++)
				{
					SWEETMath::cartesian_to_latlon(
							pos_x_D, pos_y_D, pos_z_D,
							o_pos_lon_D, o_pos_lat_D
						);

					ScalarDataArray u_D = sphereSampler.bilinear_scalar(i_dt_u_lon, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray v_D = sphereSampler.bilinear_scalar(i_dt_v_lat, o_pos_lon_D, o_pos_lat_D, true, simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_D(num_elements), vel_y_D(num_elements), vel_z_D(num_elements);
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							o_pos_lon_D, o_pos_lat_D,
							u_D, v_D,
							vel_x_D, vel_y_D, vel_z_D
						);

					// convert final points from Cartesian space to angular space
					ScalarDataArray new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-0.5*vel_x_D,
						-0.5*vel_y_D,
						-0.5*vel_z_D,

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere_geometry
					);


					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
										(pos_y_D-new_pos_y_D).reduce_maxAbs() +
										(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				/*
				 * Given the midpoint at pos_?_D, we compute the full time step
				 */

				ScalarDataArray dot2 = 2.0*(pos_x_D*pos_x_A + pos_y_D*pos_y_A + pos_z_D*pos_z_A);

				pos_x_D = dot2*pos_x_D - pos_x_A;
				pos_y_D = dot2*pos_y_D - pos_y_A;
				pos_z_D = dot2*pos_z_D - pos_z_A;


				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

			}
			else if (trajectory_method == E_TRAJECTORY_METHOD_SETTLS_HORTAL)
			{
				/**
				 * See SETTLS paper
				 * Hortal, M. (2002). The development and testing of a new two-time-level semi-Lagrangian scheme (SETTLS) in the ECMWF forecast model. Q. J. R. Meteorol. Soc., 2, 1671–1687.
				 *
				 * We use the SETTLS formulation also to compute the departure points.
				 */

				// Extrapolate velocities at departure points
				SphereData_Physical u_extrapol = 2.0*i_dt_u_lon - i_u_lon_prev;
				SphereData_Physical v_extrapol = 2.0*i_dt_v_lat - i_v_lat_prev;

				// Convert velocities along lon/lat to scalar data array
				ScalarDataArray vel_lon_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_u_lon);
				ScalarDataArray vel_lat_array = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_dt_v_lat);

				// Polar => Cartesian velocities
				ScalarDataArray vel_x_A(num_elements), vel_y_A(num_elements), vel_z_A(num_elements);
				SWEETMath::latlon_velocity_to_cartesian_velocity(
						pos_lon_A, pos_lat_A,
						vel_lon_array, vel_lat_array,
						vel_x_A, vel_y_A, vel_z_A
					);

				/*
				 * Setup iterations
				 */

				// Departure points for iterations
				ScalarDataArray pos_x_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_y_D(sphereDataConfig->physical_array_data_number_of_elements);
				ScalarDataArray pos_z_D(sphereDataConfig->physical_array_data_number_of_elements);

				doAdvectionOnSphere(
					pos_x_A,
					pos_y_A,
					pos_z_A,

					-vel_x_A,
					-vel_y_A,
					-vel_z_A,

					pos_x_D,
					pos_y_D,
					pos_z_D,

					semi_lagrangian_approximate_sphere_geometry
				);

				double diff = -1;
				for (int iters = 0; iters < semi_lagrangian_max_iterations; iters++)
				{
					SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

					/*
					 * WARNING: Never convert this to vort/div space!!!
					 * This creates some artificial waves
					 */
					ScalarDataArray vel_lon_extrapol_D = sphereSampler.bilinear_scalar(
							u_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points
						);
					ScalarDataArray vel_lat_extrapol_D = sphereSampler.bilinear_scalar(
							v_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							simVars.disc.semi_lagrangian_sampler_use_pole_pseudo_points
						);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_extrapol_D(num_elements), vel_y_extrapol_D(num_elements), vel_z_extrapol_D(num_elements);

					// polar => Cartesian coordinates
					SWEETMath::latlon_velocity_to_cartesian_velocity(
							o_pos_lon_D, o_pos_lat_D,
							vel_lon_extrapol_D, vel_lat_extrapol_D,
							vel_x_extrapol_D, vel_y_extrapol_D, vel_z_extrapol_D
						);

					ScalarDataArray new_pos_x_D(num_elements), new_pos_y_D(num_elements), new_pos_z_D(num_elements);

					doAdvectionOnSphere(
						pos_x_A,
						pos_y_A,
						pos_z_A,

						-0.5*(vel_x_extrapol_D + vel_x_A),
						-0.5*(vel_y_extrapol_D + vel_y_A),
						-0.5*(vel_z_extrapol_D + vel_z_A),

						new_pos_x_D,
						new_pos_y_D,
						new_pos_z_D,

						semi_lagrangian_approximate_sphere_geometry
					);

					if (semi_lagrangian_convergence_threshold > 0)
					{
						diff =  (pos_x_D-new_pos_x_D).reduce_maxAbs() +
								(pos_y_D-new_pos_y_D).reduce_maxAbs() +
								(pos_z_D-new_pos_z_D).reduce_maxAbs();

						if (diff < semi_lagrangian_convergence_threshold)
						{
							pos_x_D = new_pos_x_D;
							pos_y_D = new_pos_y_D;
							pos_z_D = new_pos_z_D;

							break;
						}
					}

					pos_x_D = new_pos_x_D;
					pos_y_D = new_pos_y_D;
					pos_z_D = new_pos_z_D;
				}

				if (semi_lagrangian_convergence_threshold > 0)
				{
					if (diff > semi_lagrangian_convergence_threshold)
					{
						std::cout << "WARNING: Over convergence tolerance" << std::endl;
						std::cout << "+ maxAbs: " << diff << std::endl;
						std::cout << "+ Convergence tolerance: " << semi_lagrangian_convergence_threshold << std::endl;
					}
				}

				// convert final points from Cartesian space to angular space
				SWEETMath::cartesian_to_latlon(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else
			{
				FatalError("Unknown departure point calculation method");
			}

			return;
		}

		FatalError("Only 1st and 2nd order time integration supported");
	}

};

#endif /* SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_ */
