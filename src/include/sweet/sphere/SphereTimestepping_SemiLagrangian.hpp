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

#include <sweet/ScalarDataArray.hpp>
#include <sweet/sphere/SphereData_Physical.hpp>
#include <sweet/sphere/Convert_ScalarDataArray_to_SphereDataPhysical.hpp>

#include <sweet/sphere/Convert_SphereDataPhysical_to_ScalarDataArray.hpp>
#include <sweet/sphere/SphereOperators_Sampler_SphereDataPhysical.hpp>



class SphereTimestepping_SemiLagrangian
{
	SphereOperators_Sampler_SphereDataPhysical sample2D;
	const SphereData_Config *sphereDataConfig;


public:
	SphereTimestepping_SemiLagrangian()	:
		sphereDataConfig(nullptr)
	{
	}


	void setup(
		const SphereData_Config *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sample2D.setup(sphereDataConfig);
	}


	inline
	static
	void angleToCartCoord(
			const ScalarDataArray &i_lon,
			const ScalarDataArray &i_lat,
			ScalarDataArray &o_x,
			ScalarDataArray &o_y,
			ScalarDataArray &o_z
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.number_of_elements; i++)
		{
			o_x.scalar_data[i] = std::cos(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_y.scalar_data[i] = std::sin(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_z.scalar_data[i] = std::sin(i_lat.scalar_data[i]);
		}
	}



	inline
	static
	void angleSpeedToCartVector(
			const ScalarDataArray &i_lon,
			const ScalarDataArray &i_lat,
			const ScalarDataArray &i_vel_lon,
			const ScalarDataArray &i_vel_lat,
			ScalarDataArray *o_v_x,
			ScalarDataArray *o_v_y,
			ScalarDataArray *o_v_z
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_lon.number_of_elements; i++)
		{
			o_v_x->scalar_data[i] = -i_vel_lon.scalar_data[i]*std::sin(i_lon.scalar_data[i]) - i_vel_lat.scalar_data[i]*std::cos(i_lon.scalar_data[i])*std::sin(i_lat.scalar_data[i]);
			o_v_y->scalar_data[i] = i_vel_lon.scalar_data[i]*std::cos(i_lon.scalar_data[i]) - i_vel_lat.scalar_data[i]*std::sin(i_lon.scalar_data[i])*std::sin(i_lat.scalar_data[i]);
			o_v_z->scalar_data[i] = i_vel_lat.scalar_data[i]*std::cos(i_lat.scalar_data[i]);
		}
	}



	inline
	static
	void cartToAngleCoord(
			const ScalarDataArray &i_x,
			const ScalarDataArray &i_y,
			const ScalarDataArray &i_z,
			ScalarDataArray &o_lon,
			ScalarDataArray &o_lat
	)
	{
		SWEET_THREADING_SPACE_PARALLEL_FOR
		for (std::size_t i = 0; i < i_x.number_of_elements; i++)
		{
			/*
			 * Make sure that coordinates are in valid range
			 */
#if 0
			i_x.scalar_data[i] = std::min(1., i_x.scalar_data[i]);
			i_x.scalar_data[i] = std::max(-1., i_x.scalar_data[i]);

			i_y.scalar_data[i] = std::min(1., i_y.scalar_data[i]);
			i_y.scalar_data[i] = std::max(-1., i_y.scalar_data[i]);

			i_z.scalar_data[i] = std::min(1., i_z.scalar_data[i]);
			i_z.scalar_data[i] = std::max(-1., i_z.scalar_data[i]);
#endif

			/*
			 * Now compute the angles
			 */
			o_lon.scalar_data[i] = std::atan(i_y.scalar_data[i]/i_x.scalar_data[i]);

#if SWEET_DEBUG
			if (	(std::isnan(o_lon.scalar_data[i]) != 0) ||
					(std::abs(std::isinf(o_lon.scalar_data[i])) == 1)
			)
			{
				std::cout << "Found nan/inf at position " << i << std::endl;
				std::cout << "o_lon value: " << o_lon.scalar_data[i] << std::endl;
				std::cout << "   atan(" << i_y.scalar_data[i] << ", " << i_x.scalar_data[i] << ")" << std::endl;
				FatalError("EXIT");
			}
#endif

			if (i_x.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI;
			else if (i_y.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI*2.0;

			o_lat.scalar_data[i] = std::acos(-i_z.scalar_data[i]) - M_PI*0.5;

#if SWEET_DEBUG
			if (
					(std::isnan(o_lat.scalar_data[i]) != 0) ||
					(std::abs(std::isinf(o_lat.scalar_data[i])) == 1)
			)
			{
				std::cout << "Found nan/inf at position " << i << std::endl;
				std::cout << "o_lat value: " << o_lat.scalar_data[i] << std::endl;
				std::cout << "   acos(" << i_z.scalar_data[i] << ")" << std::endl;
				FatalError("EXIT");
			}
#endif
		}
	}



	static
	double& alpha()
	{
		static double alpha = 0;
		return alpha;
	}

	double u_analytical(double i_lambda, double i_theta)
	{
		double a = 6.37122e6;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		assert(i_lambda >= 0);
		assert(i_lambda <= 2.0*M_PI);
		assert(i_theta >= -0.5*M_PI);
		assert(i_theta <= 0.5*M_PI);

		return u0*(std::cos(i_theta)*std::cos(alpha()) + std::sin(i_theta)*std::cos(i_lambda)*std::sin(alpha()));
	}


	ScalarDataArray u_analytical(
			const ScalarDataArray &i_lambda,
			const ScalarDataArray &i_theta
	)
	{
		ScalarDataArray ret;
		ret.setup(i_lambda.number_of_elements);
		for (int i = 0; i < (int)i_lambda.number_of_elements; i++)
			ret.scalar_data[i] = u_analytical(i_lambda.scalar_data[i], i_theta.scalar_data[i]);

		return ret;
	}

	double v_analytical(double i_lambda, double i_theta)
	{
		double a = 6.37122e6;
		double u0 = (2.0*M_PI*a)/(12.0*24.0*60.0*60.0);

		assert(i_lambda >= 0);
		assert(i_lambda <= 2.0*M_PI);
		assert(i_theta >= -0.5*M_PI);
		assert(i_theta <= 0.5*M_PI);

		return -u0*std::sin(i_lambda)*std::sin(alpha());
	}


	ScalarDataArray v_analytical(
			const ScalarDataArray &i_lambda,
			const ScalarDataArray &i_theta
	)
	{
		ScalarDataArray ret;
		ret.setup(i_lambda.number_of_elements);
		for (int i = 0; i < (int)i_lambda.number_of_elements; i++)
			ret.scalar_data[i] = v_analytical(i_lambda.scalar_data[i], i_theta.scalar_data[i]);

		return ret;
	}


	/**
	 * Do 1st order accurate advection on the sphere for given
	 *  - starting point,
	 *  - velocity and
	 *  - timestep size
	 * in Cartesian space 
	 */
	inline static
	void doAdvectionOnSphere(
		const ScalarDataArray &i_pos_x,
		const ScalarDataArray &i_pos_y,
		const ScalarDataArray &i_pos_z,

		const ScalarDataArray &i_vel_x,
		const ScalarDataArray &i_vel_y,
		const ScalarDataArray &i_vel_z,

		double i_dt,
		double i_inv_earth_radius,

		ScalarDataArray &o_pos_x,
		ScalarDataArray &o_pos_y,
		ScalarDataArray &o_pos_z
	)
	{
		o_pos_x = i_pos_x + i_dt*(i_vel_x)*i_inv_earth_radius;
		o_pos_y = i_pos_y + i_dt*(i_vel_y)*i_inv_earth_radius;
		o_pos_z = i_pos_z + i_dt*(i_vel_z)*i_inv_earth_radius;

		ScalarDataArray norm = (o_pos_x*o_pos_x + o_pos_y*o_pos_y + o_pos_z*o_pos_z).inv_sqrt();

		o_pos_x *= norm;
		o_pos_y *= norm;
		o_pos_z *= norm;
	}



	void semi_lag_departure_points_settls(
			const SphereData_Physical &i_u_lon_prev,	// Velocities at time t-1
			const SphereData_Physical &i_v_lat_prev,

			const SphereData_Physical &i_u_lon, 		// Velocities at time t
			const SphereData_Physical &i_v_lat,

			const ScalarDataArray &i_pos_lon_a,		// Position of arrival points lon/lat
			const ScalarDataArray &i_pos_lat_a,

			double i_dt,				///< time step size
			double i_earth_radius,

			ScalarDataArray &o_pos_lon_d, 	///< Position of departure points x / y
			ScalarDataArray &o_pos_lat_d,

			int i_timestepping_order,
			int max_iters = 10,
			double i_convergence_tolerance = 1e-8,
			int i_approximate_sphere_geometry = 0,
			bool use_interpolation_limiters = false
	)
	{
#if 0
#if 1
		i_u_lon_prev.file_read_raw("i_u_lon_prev.raw");
		i_v_lat_prev.file_read_raw("i_v_lat_prev.raw");

		i_u_lon.file_read_raw("i_u_lon.raw");
		i_v_lat.file_read_raw("i_v_lat.raw");

//		i_pos_lon_a.file_read_raw("i_pos_lon_a.raw");
//		i_pos_lat_a.file_read_raw("i_pos_lat_a.raw");
#else
		i_u_lon_prev.file_write_raw("i_u_lon_prev.raw");
		i_v_lat_prev.file_write_raw("i_v_lat_prev.raw");

		i_u_lon.file_write_raw("i_u_lon.raw");
		i_v_lat.file_write_raw("i_v_lat.raw");

		i_pos_lon_a.file_write_raw("i_pos_lon_a.raw");
		i_pos_lat_a.file_write_raw("i_pos_lat_a.raw");
#endif
#endif


#if SWEET_DEBUG
		if (i_u_lon_prev.physical_isAnyNaNorInf())
			FatalError("start NaN/Inf in i_u_lon_prev");

		if (i_v_lat_prev.physical_isAnyNaNorInf())
			FatalError("start NaN/Inf in i_v_lat_prev");

		if (i_u_lon.physical_isAnyNaNorInf())
			FatalError("start NaN/Inf in i_u_lon");

		if (i_v_lat.physical_isAnyNaNorInf())
			FatalError("start NaN/Inf in i_v_lat");
#endif

		if (i_approximate_sphere_geometry == 0)
		{
			FatalError("TODO: Implement me: i_approximate_sphere_geometry==0 to avoid approximation of SL on sphere");
		}

		std::size_t num_elements = i_pos_lon_a.number_of_elements;
		double inv_earth_radius = 1.0/i_earth_radius;

		if (i_timestepping_order == 1)
		{
			ScalarDataArray pos_x_a(num_elements);
			ScalarDataArray pos_y_a(num_elements);
			ScalarDataArray pos_z_a(num_elements);

			// polar => Cartesian coordinates
			angleToCartCoord(
					i_pos_lon_a, i_pos_lat_a,
					pos_x_a, pos_y_a, pos_z_a
				);

			ScalarDataArray u_lon = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_u_lon);
			ScalarDataArray v_lat = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_v_lat);

			ScalarDataArray vel_x(num_elements);
			ScalarDataArray vel_y(num_elements);
			ScalarDataArray vel_z(num_elements);

			// polar => Cartesian coordinates
			angleSpeedToCartVector(
					i_pos_lon_a, i_pos_lat_a,
					u_lon, v_lat,
					&vel_x, &vel_y, &vel_z
				);

			// go to departure point
			ScalarDataArray pos_x_d = pos_x_a - vel_x*i_dt*inv_earth_radius;
			ScalarDataArray pos_y_d = pos_y_a - vel_y*i_dt*inv_earth_radius;
			ScalarDataArray pos_z_d = pos_z_a - vel_z*i_dt*inv_earth_radius;

			// normalize
			ScalarDataArray norm = (pos_x_d*pos_x_d + pos_y_d*pos_y_d + pos_z_d*pos_z_d).inv_sqrt();

			pos_x_d *= norm;
			pos_y_d *= norm;
			pos_z_d *= norm;

			cartToAngleCoord(pos_x_d, pos_y_d, pos_z_d, o_pos_lon_d, o_pos_lat_d);
			return;
		}

		if (i_timestepping_order == 2)
		{
			/**
			 * See SETTLS paper
			 * Hortal, M. (2002). The development and testing of a new two-time-level semi-Lagrangian scheme (SETTLS) in the ECMWF forecast model. Q. J. R. Meteorol. Soc., 2, 1671–1687.
			 * 
			 * Note, that the equations are rearranged, see 
			 */

			// Extrapolate velocities at departure points
			SphereData_Physical u_extrapol = 2.0*i_u_lon - i_u_lon_prev;
			SphereData_Physical v_extrapol = 2.0*i_v_lat - i_v_lat_prev;

			// Compute Cartesian arrival points
			ScalarDataArray pos_x_a(num_elements);
			ScalarDataArray pos_y_a(num_elements);
			ScalarDataArray pos_z_a(num_elements);
			angleToCartCoord(
					i_pos_lon_a, i_pos_lat_a,
					pos_x_a, pos_y_a, pos_z_a
				);

			// Convert velocities along lon/lat to scalardata array
			ScalarDataArray u_lon = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_u_lon);
			ScalarDataArray v_lat = Convert_SphereDataPhysical_to_ScalarDataArray::physical_convert(i_v_lat);

			// Compute Cartesian velocities
			ScalarDataArray vel_x(num_elements);
			ScalarDataArray vel_y(num_elements);
			ScalarDataArray vel_z(num_elements);

			// Polar => Cartesian coordinates
			angleSpeedToCartVector(
					i_pos_lon_a, i_pos_lat_a,
					u_lon, v_lat,
					&vel_x, &vel_y, &vel_z
				);

			/*
			 * Setup iterations
			 */
			// Departure points for iterations
			ScalarDataArray pos_x_d = pos_x_a;
			ScalarDataArray pos_y_d = pos_y_a;
			ScalarDataArray pos_z_d = pos_z_a;


			double diff = 999;
			int iters = 0;
			for (; iters < max_iters; iters++)
			{
				cartToAngleCoord(pos_x_d, pos_y_d, pos_z_d, o_pos_lon_d, o_pos_lat_d);

#if SWEET_DEBUG
				bool stop = false;
				if (pos_x_d.reduce_isAnyNaNorInf())
				{
					std::cout << "iters: " << iters << std::endl;
					std::cout << "NaN/Inf in pos_x_d" << std::endl;
					stop = true;
				}

				if (pos_y_d.reduce_isAnyNaNorInf())
				{
					std::cout << "iters: " << iters << std::endl;
					std::cout << "NaN/Inf in pos_y_d" << std::endl;
					stop = true;
				}

				if (pos_z_d.reduce_isAnyNaNorInf())
				{
					std::cout << "iters: " << iters << std::endl;
					std::cout << "NaN/Inf in pos_z_d" << std::endl;
					stop = true;
				}


				if (o_pos_lon_d.reduce_isAnyNaNorInf())
				{
					std::cout << "iters: " << iters << std::endl;
					std::cout << "NaN/Inf in o_pos_lon_d" << std::endl;
					stop = true;
				}

				if (o_pos_lat_d.reduce_isAnyNaNorInf())
				{
					std::cout << "iters: " << iters << std::endl;
					std::cout << "NaN/Inf in o_pos_lat_d" << std::endl;
					stop = true;
				}

				if (stop)
				{
					FatalError("STOP");
				}
#endif


#if 0
				// TODO: This should be linear
				// TODO: This should be linear
				// TODO: This should be linear

				ScalarDataArray u_lon_extrapol = sample2D.bicubic_scalar(u_extrapol, o_pos_lon_d, o_pos_lat_d, true, use_interpolation_limiters);
				ScalarDataArray v_lat_extrapol = sample2D.bicubic_scalar(v_extrapol, o_pos_lon_d, o_pos_lat_d, true, use_interpolation_limiters);
#else
				ScalarDataArray u_lon_extrapol = sample2D.bilinear_scalar(u_extrapol, o_pos_lon_d, o_pos_lat_d, true);
				ScalarDataArray v_lat_extrapol = sample2D.bilinear_scalar(v_extrapol, o_pos_lon_d, o_pos_lat_d, true);
#endif

				// convert extrapolated velocities to Cartesian velocities
				ScalarDataArray vel_x_extrapol(num_elements);
				ScalarDataArray vel_y_extrapol(num_elements);
				ScalarDataArray vel_z_extrapol(num_elements);

				// polar => Cartesian coordinates
				angleSpeedToCartVector(
						o_pos_lon_d, o_pos_lat_d,
						u_lon_extrapol, v_lat_extrapol,
						&vel_x_extrapol, &vel_y_extrapol, &vel_z_extrapol
					);

				// pos_x_d = pos_x_a - i_dt*0.5*(vel_x_extrapol + vel_x)*inv_earth_radius;
				// pos_y_d = pos_y_a - i_dt*0.5*(vel_y_extrapol + vel_y)*inv_earth_radius;
				// pos_z_d = pos_z_a - i_dt*0.5*(vel_z_extrapol + vel_z)*inv_earth_radius;
				ScalarDataArray new_pos_x_d(num_elements);
				ScalarDataArray new_pos_y_d(num_elements);
				ScalarDataArray new_pos_z_d(num_elements);
				
				doAdvectionOnSphere(
					pos_x_a,
					pos_y_a,
					pos_z_a,

					vel_x_extrapol + vel_x,
					vel_y_extrapol + vel_y,
					vel_z_extrapol + vel_z,

					-i_dt*0.5,
					inv_earth_radius,

					new_pos_x_d,
					new_pos_y_d,
					new_pos_z_d
				);

				diff =  (pos_x_d-new_pos_x_d).reduce_maxAbs() +
						(pos_y_d-new_pos_y_d).reduce_maxAbs() +
						(pos_z_d-new_pos_z_d).reduce_maxAbs();

				pos_x_d = new_pos_x_d;
				pos_y_d = new_pos_y_d;
				pos_z_d = new_pos_z_d;

				if (diff < i_convergence_tolerance)
				   break;
			}

			if (diff > i_convergence_tolerance)
			{
				std::cout << "WARNING: Over convergence tolerance" << std::endl;
				std::cout << "+ Iterations: " << iters << std::endl;
				std::cout << "+ maxAbs: " << diff << std::endl;
				std::cout << "+ Convergence tolerance: " << i_convergence_tolerance << std::endl;
			}

			// convert final points from Cartesian space to angular space
			cartToAngleCoord(pos_x_d, pos_y_d, pos_z_d, o_pos_lon_d, o_pos_lat_d);
			return;
		}

		FatalError("Only 1st and 2nd order time integration supported");
	}

};

#endif /* SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_ */
