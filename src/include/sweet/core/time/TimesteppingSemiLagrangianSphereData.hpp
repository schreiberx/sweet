/*
 * SphereDataSemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 *  Updated to sphere on 28th March 2018
 */
#ifndef SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_
#define SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_

#include <limits>
#include <sweet/core/ScalarDataArray.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/Convert_ScalarDataArray_to_SphereDataPhysical.hpp>
#include <sweet/core/sphere/Convert_SphereDataPhysical_to_ScalarDataArray.hpp>
#include <sweet/core/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>
#include <sweet/core/VectorMath.hpp>
#include "ShackTimesteppingSemiLagrangianSphereData.hpp"

namespace sweet
{

class SphereTimestepping_SemiLagrangian
{
	const SphereDataConfig *sphereDataConfig;
	sweet::ShackTimesteppingSemiLagrangianSphereData *slData;

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
			ShackTimesteppingSemiLagrangianSphereData *i_slData
	)	:
		sphereDataConfig(nullptr),
		slData(i_slData)
	{
	}


	SphereTimestepping_SemiLagrangian()	:
		sphereDataConfig(nullptr),
		slData(nullptr)
	{
	}


	void setup(
			const SphereDataConfig *i_sphereDataConfig,
			ShackTimesteppingSemiLagrangianSphereData *i_slData,
			int i_timestepping_order
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sphereSampler.setup(sphereDataConfig);

		if (slData->semi_lagrangian_departure_point_method == "settls")
		{
			trajectory_method = E_TRAJECTORY_METHOD_SETTLS_HORTAL;
		}
		else if (slData->semi_lagrangian_departure_point_method == "canonical")
		{
			trajectory_method = E_TRAJECTORY_METHOD_CANONICAL;
		}
		else if (slData->semi_lagrangian_departure_point_method == "midpoint")
		{
			trajectory_method = E_TRAJECTORY_METHOD_MIDPOINT_RITCHIE;
		}
		else
		{
			SWEETError(std::string("Trajectory method '")+slData->semi_lagrangian_departure_point_method+"' not supported");
		}


		timestepping_order = i_timestepping_order;
		semi_lagrangian_max_iterations = slData->semi_lagrangian_max_iterations;
		semi_lagrangian_convergence_threshold = slData->semi_lagrangian_convergence_threshold;
		semi_lagrangian_approximate_sphere_geometry = slData->semi_lagrangian_approximate_sphere_geometry;
		semi_lagrangian_interpolation_limiter = slData->semi_lagrangian_interpolation_limiter;



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

		VectorMath::point_latlon_to_cartesian__array(
				pos_lon_A, pos_lat_A,
				pos_x_A, pos_y_A, pos_z_A
			);


#if 0
		//SphereData_Physical u_lon_
		if (simVars.sim.sphere_use_fsphere)
		{
			error.set("Not supported");
		}

		sl_coriolis.setup_if_required(sphereDataConfig);
		sl_coriolis.physical_update_lambda_gaussian_grid(
			[&](double lon, double mu, double &o_data)
			{
				o_data = mu*2.0*simVars.sim.sphere_rotating_coriolis_omega*simVars.sim.sphere_radius;
			}
		);
#endif

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
			std::cout << "i_approximate_sphere_geometry" << std::endl;
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
			VectorMath::normalize(o_pos_x, o_pos_y, o_pos_z);
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

			VectorMath::cross_prod(
					i_pos_x, i_pos_y, i_pos_z,
					i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z,
					rotation_axis_x, rotation_axis_y, rotation_axis_z
				);
#if 1
			/*
			 * Normalize rotation axis since it's likely not normalized yet
			 */

			static int asdf = 1;

			if (asdf == 1)
			{
				std::cout << "TODO: REPLACE ME" << std::endl;
				std::cout << "TODO: REPLACE ME" << std::endl;
				std::cout << "TODO: REPLACE ME" << std::endl;
				std::cout << "TODO: REPLACE ME" << std::endl;
				std::cout << "TODO: REPLACE ME" << std::endl;
				std::cout << "TODO: REPLACE ME" << std::endl;
				asdf = 0;
			}

			/*
			 * TODO: replace this!
			 *
			 * Use a formulation without normalization by using the
			 * angular vector / velocity e.g. by using quaternions
			 */
			VectorMath::normalize_with_threshold(
					rotation_axis_x,
					rotation_axis_y,
					rotation_axis_z
				);

			/*
			 * Compute rotation angle
			 *
			 * No rescaling by 1/(2pi) since the change of angle is directly
			 * given by the magnitude of the angular velocity by its definition
			 */
			ScalarDataArray angle = VectorMath::length(i_dt_velocity_x, i_dt_velocity_y, i_dt_velocity_z);
#else
			// doesn't work
			ScalarDataArray angle = i_dt_velocity_x;
			angle.physical_set_all(1.0);
#endif
			/*
			 * Rotate
			 */
			VectorMath::point_rotate_3d_normalized_rotation_axis__array(
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



	/*
	 * Compute SL departure points on unit sphere for given dt*(u,v) velocities
	 *
	 * All this is for the unit sphere and unit time!
	 *
	 * Hence, the velocities need to be rescaled by "dt/radius"
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
			VectorMath::velocity_latlon_to_cartesian__array(
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
			VectorMath::point_cartesian_to_latlon__array(
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
				VectorMath::velocity_latlon_to_cartesian__array(
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

					VectorMath::point_cartesian_to_latlon__array(
							pos_x_mid, pos_y_mid, pos_z_mid,
							pos_lon_mid, pos_lat_mid
						);

					ScalarDataArray vel_u_mid = sphereSampler.bilinear_scalar(vel_lon, pos_lon_mid, pos_lat_mid, true, slData->semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray vel_v_mid = sphereSampler.bilinear_scalar(vel_lat, pos_lon_mid, pos_lat_mid, true, slData->semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_mid(num_elements), vel_y_mid(num_elements), vel_z_mid(num_elements);
					VectorMath::velocity_latlon_to_cartesian__array(
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
				VectorMath::point_cartesian_to_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
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
				VectorMath::velocity_latlon_to_cartesian__array(
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
					VectorMath::point_cartesian_to_latlon__array(
							pos_x_D, pos_y_D, pos_z_D,
							o_pos_lon_D, o_pos_lat_D
						);

					ScalarDataArray u_D = sphereSampler.bilinear_scalar(i_dt_u_lon, o_pos_lon_D, o_pos_lat_D, true, slData->semi_lagrangian_sampler_use_pole_pseudo_points);
					ScalarDataArray v_D = sphereSampler.bilinear_scalar(i_dt_v_lat, o_pos_lon_D, o_pos_lat_D, true, slData->semi_lagrangian_sampler_use_pole_pseudo_points);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_D(num_elements), vel_y_D(num_elements), vel_z_D(num_elements);
					VectorMath::velocity_latlon_to_cartesian__array(
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


				VectorMath::point_cartesian_to_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

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
				VectorMath::velocity_latlon_to_cartesian__array(
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
					VectorMath::point_cartesian_to_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);

					/*
					 * WARNING: Never convert this to vort/div space!!!
					 * This creates some artificial waves
					 */
					ScalarDataArray vel_lon_extrapol_D = sphereSampler.bilinear_scalar(
							u_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							slData->semi_lagrangian_sampler_use_pole_pseudo_points
						);
					ScalarDataArray vel_lat_extrapol_D = sphereSampler.bilinear_scalar(
							v_extrapol,
							o_pos_lon_D, o_pos_lat_D,
							true,
							slData->semi_lagrangian_sampler_use_pole_pseudo_points
						);

					// convert extrapolated velocities to Cartesian velocities
					ScalarDataArray vel_x_extrapol_D(num_elements), vel_y_extrapol_D(num_elements), vel_z_extrapol_D(num_elements);

					// polar => Cartesian coordinates
					VectorMath::velocity_latlon_to_cartesian__array(
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
				VectorMath::point_cartesian_to_latlon__array(pos_x_D, pos_y_D, pos_z_D, o_pos_lon_D, o_pos_lat_D);
			}
			else
			{
				SWEETError("Unknown departure point calculation method");
			}

			return;
		}

		SWEETError("Only 1st and 2nd order time integration supported");
	}



	/*
	 * Interpolation of prognostic fields at departure points.
	 *
	 * We assume the velocity U-V to be the SL advected field!
	 */
	void apply_sl_timeintegration_vd(
			const SphereOperators &i_ops,

			const SphereData_Spectral &i_phi,
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,

			const ScalarDataArray &i_pos_lon_d,
			const ScalarDataArray &i_pos_lat_d,

			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	)
	{
		o_phi.setup_if_required(i_phi.sphereDataConfig);
		o_vrt.setup_if_required(i_phi.sphereDataConfig);
		o_div.setup_if_required(i_phi.sphereDataConfig);

		o_phi = sphereSampler.bicubic_scalar_ret_phys(
				i_phi.toPhys(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);

		o_vrt = sphereSampler.bicubic_scalar_ret_phys(
				i_vrt.toPhys(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);

		o_div = sphereSampler.bicubic_scalar_ret_phys(
				i_div.toPhys(),
				i_pos_lon_d, i_pos_lat_d,
				false,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);
	}


	/*
	 * Interpolation of prognostic fields at departure points.
	 *
	 * We assume the velocity U-V to be the SL advected field!
	 */
	void apply_sl_timeintegration_uv(
			const SphereOperators &i_ops,

			const SphereData_Spectral &i_phi,
			const SphereData_Spectral &i_vrt,
			const SphereData_Spectral &i_div,

			const ScalarDataArray &i_pos_lon_D,
			const ScalarDataArray &i_pos_lat_D,

			SphereData_Spectral &o_phi,
			SphereData_Spectral &o_vrt,
			SphereData_Spectral &o_div
	)
	{
		o_phi.setup_if_required(i_phi.sphereDataConfig);
		o_vrt.setup_if_required(i_phi.sphereDataConfig);
		o_div.setup_if_required(i_phi.sphereDataConfig);


		/*************************************************************************
		 * Phi
		 *************************************************************************
		 */
		o_phi = sphereSampler.bicubic_scalar_ret_phys(
				i_phi.toPhys(),
				i_pos_lon_D, i_pos_lat_D,
				false,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);


	#if 1

		/*************************************************************************
		 * Prepare rotation system for handling velocities
		 *************************************************************************
		 */

		ScalarDataArray P_x_D, P_y_D, P_z_D;
		VectorMath::point_latlon_to_cartesian__array(i_pos_lon_D, i_pos_lat_D, P_x_D, P_y_D, P_z_D);

		ScalarDataArray	&P_x_A = pos_x_A,
						&P_y_A = pos_y_A,
						&P_z_A = pos_z_A;

		/*
		 * Compute rotation angle
		 */

		ScalarDataArray rotation_angle_ =
				VectorMath::dot_prod(
					P_x_D, P_y_D, P_z_D,
					P_x_A, P_y_A, P_z_A
				);

		// Can be slightly larger than 1, leading to NaN, hence this hack
		rotation_angle_ = VectorMath::min(rotation_angle_, 1.0);

		ScalarDataArray rotation_angle = VectorMath::arccos(rotation_angle_);


		/*
		 * Compute Rotation axis
		 */
		ScalarDataArray rot_x, rot_y, rot_z;
		VectorMath::cross_prod(
				P_x_D, P_y_D, P_z_D,
				P_x_A, P_y_A, P_z_A,
				rot_x, rot_y, rot_z
			);

		VectorMath::normalize_with_threshold(rot_x, rot_y, rot_z);



		/*************************************************************************
		 * Velocity
		 *************************************************************************
		 */

		SphereData_Physical u_tmp, v_tmp;
		i_ops.vrtdiv_to_uv(i_vrt, i_div, u_tmp, v_tmp);

		SphereData_Physical u_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
				u_tmp,
				i_pos_lon_D, i_pos_lat_D,
				true,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);

		SphereData_Physical v_tmp_D = sphereSampler.bicubic_scalar_ret_phys(
				v_tmp,
				i_pos_lon_D, i_pos_lat_D,
				true,
				slData->semi_lagrangian_sampler_use_pole_pseudo_points,
				slData->semi_lagrangian_interpolation_limiter
			);

		/*
		 * Convert to Cartesian space
		 */

		ScalarDataArray V_lon_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(u_tmp_D);
		ScalarDataArray V_lat_D = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(v_tmp_D);

		ScalarDataArray V_x_D, V_y_D, V_z_D;
		VectorMath::velocity_latlon_to_cartesian__array(
				i_pos_lon_D,
				i_pos_lat_D,
				V_lon_D,
				V_lat_D,
				V_x_D,
				V_y_D,
				V_z_D
			);

		/*
		 * Rotate to velocity vector
		 */
		ScalarDataArray V_x_A, V_y_A, V_z_A;
		VectorMath::vector_rotate_3d_normalized_rotation_axis__array(
				V_x_D, V_y_D, V_z_D,
				rotation_angle,
				rot_x, rot_y, rot_z,
				V_x_A, V_y_A, V_z_A
			);


		/*
		 * Return velocity in lat/lon space
		 */
		ScalarDataArray V_lon_A, V_lat_A;
		VectorMath::velocity_cartesian_to_latlon__array(
				pos_lon_A,
				pos_lat_A,
				V_x_A,
				V_y_A,
				V_z_A,
				V_lon_A, V_lat_A
		);

		i_ops.uv_to_vrtdiv(
				Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lon_A, i_vrt.sphereDataConfig),
				Convert_ScalarDataArray_to_SphereDataPhysical::convert(V_lat_A, i_vrt.sphereDataConfig),
				o_vrt, o_div
			);

	#else

		op.uv_to_vrtdiv(u_tmp_D, v_tmp_D, o_vrt, o_div);

	#endif
	}

};

}

#endif
