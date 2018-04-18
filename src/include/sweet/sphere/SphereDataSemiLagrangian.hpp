/*
 * SphereDataSemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 *
 *  Updated to sphere on 28th March 2018
 */
#ifndef SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_
#define SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_

#include <sweet/sphere/Convert_SphereData_to_ScalarDataArray.hpp>
#include <sweet/sphere/Convert_SphereDataPhysical_to_ScalarDataArray.hpp>
#include <sweet/sphere/Convert_ScalarDataArray_to_SphereData.hpp>
#include <sweet/sphere/SphereStaggering.hpp>
#include <sweet/sphere/SphereDataSampler.hpp>
#include "SphereData.hpp"
#include <sweet/ScalarDataArray.hpp>



class SphereDataSemiLagrangian
{
	SphereDataSampler sample2D;
	const SphereDataConfig *sphereDataConfig;


public:
	SphereDataSemiLagrangian()	:
		sphereDataConfig(nullptr)
	{
	}


	void setup(
		double i_domain_size[2],
		const SphereDataConfig *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sample2D.setup(sphereDataConfig);
	}


#if 0
	inline
	static
	void angleToCartCoord(
			double i_lon,
			double i_lat,
			double *o_x
	)
	{
		i_lat += M_PI*0.5;
		o_x[0] = std::sin(i_lat)*std::cos(i_lon);
		o_x[1] = std::sin(i_lat)*std::sin(i_lon);
		o_x[2] = -std::cos(i_lat);
	}
#else

	inline
	static
	void angleToCartCoord(
			double i_lambda,
			double i_theta,
			double *o_x
	)
	{
		o_x[0] = std::cos(i_lambda)*std::cos(i_theta);
		o_x[1] = std::sin(i_lambda)*std::cos(i_theta);
		o_x[2] = std::sin(i_theta);
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < i_lon.number_of_elements; i++)
		{
			o_x.scalar_data[i] = std::cos(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_y.scalar_data[i] = std::sin(i_lon.scalar_data[i])*std::cos(i_lat.scalar_data[i]);
			o_z.scalar_data[i] = std::sin(i_lat.scalar_data[i]);
		}
	}

#endif

	inline
	static
	void angleSpeedToCartVector(
			double i_lambda,
			double i_theta,
			double i_u,
			double i_v,
			double *o_x
	)
	{
		//i_v = -i_v;
		o_x[0] = -i_u*std::sin(i_lambda) - i_v*std::cos(i_lambda)*std::sin(i_theta);
		o_x[1] = i_u*std::cos(i_lambda) - i_v*std::sin(i_lambda)*std::sin(i_theta);
		o_x[2] = i_v*std::cos(i_theta);
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
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
			const double i_x[3],
			double *o_lon,
			double *o_lat
	)
	{
#if 0
		*o_lon = std::acos(i_x[0]/std::sqrt(i_x[0]*i_x[0] + i_x[1]*i_x[1]));
		if (i_x[1] < 0)
			*o_lon = 2.0*M_PI-*o_lon;
#else

		*o_lon = std::atan(i_x[1]/i_x[0]);

		if (i_x[0] < 0)
			*o_lon += M_PI;
		else if (i_x[1] < 0)
			*o_lon += M_PI*2.0;

#endif

		*o_lat = std::acos(-i_x[2]);
		*o_lat -= M_PI*0.5;
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
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < i_x.number_of_elements; i++)
		{
			o_lon.scalar_data[i] = std::atan(i_y.scalar_data[i]/i_x.scalar_data[i]);

			if (i_x.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI;
			else if (i_y.scalar_data[i] < 0)
				o_lon.scalar_data[i] += M_PI*2.0;

			o_lat.scalar_data[i] = std::acos(-i_z.scalar_data[i]) - M_PI*0.5;
		}
	}


	inline
	static
	double length(
			const double *i_x
	)
	{
		return std::sqrt(i_x[0]*i_x[0] + i_x[1]*i_x[1] + i_x[2]*i_x[2]);
	}


#if 0
	inline
	static
	void angleToTangentialSpace(
			double i_lon,
			double i_lat,
			double *o_x
	)
	{
		angleToCartCoord(i_lon, i_lat, o_x);

		o_x[3+0] = -std::sin(i_lon);
		o_x[3+1] = std::cos(i_lon);
		o_x[3+2] = 0;

		cross(&o_x[0], &o_x[3], &o_x[6]);
	}


	inline
	static
	void cross(
			const double *i_x,
			const double *i_y,
			double *o_z
	)
	{
		o_z[0] = i_x[1]*i_y[2] - i_x[2]*i_y[1];
		o_z[1] = i_x[2]*i_y[0] - i_x[0]*i_y[2];
		o_z[2] = i_x[0]*i_y[1] - i_x[1]*i_y[0];
	}

#endif


#if 0
	/**
	 * Add a vector on the sphere to a position
	 */
	void sphereCoordMinusSurfaceVector(
			const ScalarDataArray &i_pos_lon,	/// longitude angle
			const ScalarDataArray &i_pos_lat,	/// latitude angle

			const ScalarDataArray &i_vec_lon,	/// velocity along longitude
			const ScalarDataArray &i_vec_lat,	/// velocity along latitude

			double i_dt,

			ScalarDataArray &o_pos_lon,
			ScalarDataArray &o_pos_lat
	)
	{
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < i_pos_lon.number_of_elements; i++)
		{
#if 1
			double lambda = i_pos_lon.scalar_data[i];
			double theta = i_pos_lat.scalar_data[i];

			double u = i_vec_lon.scalar_data[i];
			double v = i_vec_lat.scalar_data[i];

			double XG[3];
			angleToCartCoord(lambda, theta, XG);

			double XVK[3];
			angleSpeedToCartVector(lambda, theta, u, v, XVK);

			double bk = 1.0/std::sqrt(
								1.0
								+ i_dt*i_dt*(XVK[0]*XVK[0] + XVK[1]*XVK[1] + XVK[2]*XVK[2])
								- 2.0*i_dt*(XVK[0]*XG[0] + XVK[1]*XG[1] + XVK[2]*XG[2])
							);

			double XKn[3];

			XKn[0] = bk*(XG[0] - i_dt*XVK[0]);
			XKn[1] = bk*(XG[1] - i_dt*XVK[1]);
			XKn[2] = bk*(XG[2] - i_dt*XVK[2]);

			double lambda_n = std::atan(XKn[1]/XKn[0]);
			double theta_n = std::asin(XKn[2]);

			cartToAngleCoord(XKn, &lambda_n, &theta_n);

			o_pos_lon.scalar_data[i] = lambda_n;
			o_pos_lat.scalar_data[i] = theta_n;

#elif 1

			// compute 3D pos
			double tanBasis[9];
			angleToTangentialSpace(i_pos_lon.scalar_data[i], i_pos_lat.scalar_data[i], tanBasis);

			// compute 3D velocity
			double velVec[3];
			for (int j = 0; j < 3; j++)
				velVec[j] = (tanBasis[j+3]*i_vec_lon.scalar_data[i] + tanBasis[j+6]*i_vec_lat.scalar_data[i])*i_dt;

#if 0
			// compute rotating advection along longitude

			double velBasis[9];
			velBasis[0] = 0;
			velBasis[1] = 0;
			velBasis[2] = 1;

			velBasis[3] = tanBasis[0];
			velBasis[4] = tanBasis[1];
			velBasis[5] = 0;//tanBasis[2];

			velBasis[6] = tanBasis[1];
			velBasis[7] = -tanBasis[0];
			velBasis[8] = 0;//tanBasis[2];

			// Follow vel along unit circle
			// start at relative position (1,0)

			double vel_length = i_vec_lon.scalar_data[i];
			double pos[2];
			pos[0] = std::cos(vel_length)-1.0;
			pos[1] = std::sin(vel_length);

			// update velVec
			velVec[0] = 0;
			for (int i = 0; i < 3; i++)
				velVec[0] += velBasis[3+i*3]*pos[i];

			velVec[1] = 0;
			for (int i = 0; i < 3; i++)
				velVec[1] += velBasis[6+i*3]*pos[i];

			velVec[2] = 0;

			// add
			double finalPos[3];
			for (int j = 0; j < 3; j++)
				finalPos[j] = tanBasis[j] + velVec[j];

#elif 0
			/*
			 * Project vector accurately on sphere
			 *
			 * Not fully tested yet, but seems to work betterh than other approach!!!
			 */

			// velocity basis
			double velBasis[9];
			for (int j = 0; j < 3; j++)
				velBasis[j] = tanBasis[j];

			double vel_length = std::sqrt(i_vec_lon.scalar_data[i]*i_vec_lon.scalar_data[i] + i_vec_lat.scalar_data[i]*i_vec_lat.scalar_data[i])*i_dt;
			double inv_vel_len = 1.0/vel_length;
			for (int j = 0; j < 3; j++)
				velBasis[3+j] = velVec[j]*inv_vel_len;

			cross(&velBasis[0], &velBasis[3], &velBasis[6]);

			// Follow vel along unit circle
			// start at relative position (1,0)
			double pos[2];
			pos[0] = std::cos(vel_length)-1.0;
			pos[1] = std::sin(vel_length);

			// update velVec
			velVec[0] = 0;
			for (int i = 0; i < 2; i++)
				velVec[0] += velBasis[i*3]*pos[i];

			velVec[1] = 0;
			for (int i = 0; i < 2; i++)
				velVec[1] += velBasis[3+i*3]*pos[i];

			// add
			double finalPos[3];
			for (int j = 0; j < 3; j++)
				finalPos[j] = tanBasis[j] + velVec[j];

#else

			// add
			double finalPos[3];
			for (int j = 0; j < 3; j++)
				finalPos[j] = tanBasis[j] + velVec[j];

			// normalize
			double inv_len = 1.0/length(finalPos);
			for (int j = 0; j < 3; j++)
				finalPos[j] *= inv_len;

#endif

			cartToAngleCoord(finalPos, &o_pos_lon.scalar_data[i], &o_pos_lat.scalar_data[i]);

#else

			double vx = i_vec_lon.scalar_data[i]*i_dt;
			double vy = i_vec_lat.scalar_data[i]*i_dt;

			double px = i_pos_lon.scalar_data[i];
			double py = i_pos_lat.scalar_data[i];

			assert(px >= 0);
			assert(px <= 2.0*M_PI);
			assert(py >= -0.5*M_PI);
			assert(py <= 0.5*M_PI);

			o_pos_lon.scalar_data[i] = px + vx/std::cos(py);
			o_pos_lat.scalar_data[i] = py + vy;

#if 0
			if (o_pos_lat.scalar_data[i] > M_PI*0.5)
			{
				o_pos_lat.scalar_data[i] = M_PI - o_pos_lat.scalar_data[i];
				o_pos_lon.scalar_data[i] += M_PI;
				if (o_pos_lon.scalar_data[i] > 2.0*M_PI)
					o_pos_lon.scalar_data[i] -= 2.0*M_PI;
			}

			if (o_pos_lat.scalar_data[i] < -M_PI*0.5)
			{
				o_pos_lat.scalar_data[i] = -M_PI - o_pos_lat.scalar_data[i];
				o_pos_lon.scalar_data[i] += M_PI;
				if (o_pos_lon.scalar_data[i] > 2.0*M_PI)
					o_pos_lon.scalar_data[i] -= 2.0*M_PI;
			}
#endif

#endif
		}
	}
#endif



#if 1
	//double alpha = 1.5708;
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
#endif



#if 1
	void semi_lag_departure_points_settls(
			const SphereDataPhysical &i_u_lon_prev,	// Velocities at time t-1
			const SphereDataPhysical &i_v_lat_prev,

			const SphereDataPhysical &i_u_lon, 		// Velocities at time t
			const SphereDataPhysical &i_v_lat,

			const ScalarDataArray &i_pos_lon_a,	// Position of arrival points lon/lat
			const ScalarDataArray &i_pos_lat_a,

			double i_dt,				///< time step size
			double i_earth_radius,

			ScalarDataArray &o_pos_lon_d, 	///< Position of departure points x / y
			ScalarDataArray &o_pos_lat_d,

			int i_timestepping_order,
			int max_iters = 2,
			double i_convergence_tolerance = 1e-5
	)
	{
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

			ScalarDataArray u_lon = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_u_lon);
			ScalarDataArray v_lat = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_v_lat);

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
			// Extrapolate velocities at departure points
			SphereDataPhysical u_extrapol = 2.0*i_u_lon - i_u_lon_prev;
			SphereDataPhysical v_extrapol = 2.0*i_v_lat - i_v_lat_prev;

			// Compute cartesian arrival points
			ScalarDataArray pos_x_a(num_elements);
			ScalarDataArray pos_y_a(num_elements);
			ScalarDataArray pos_z_a(num_elements);
			angleToCartCoord(
					i_pos_lon_a, i_pos_lat_a,
					pos_x_a, pos_y_a, pos_z_a
				);

			// convert velocities along lon/lat to scalardata array
			ScalarDataArray u_lon = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_u_lon);
			ScalarDataArray v_lat = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_v_lat);

			// compute Cartesian velocities
			ScalarDataArray vel_x(num_elements);
			ScalarDataArray vel_y(num_elements);
			ScalarDataArray vel_z(num_elements);

			// polar => Cartesian coordinates
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

				ScalarDataArray u_lon_extrapol = sample2D.bilinear_scalar(u_extrapol, o_pos_lon_d, o_pos_lat_d, true);
				ScalarDataArray v_lat_extrapol = sample2D.bilinear_scalar(v_extrapol, o_pos_lon_d, o_pos_lat_d, true);

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

				pos_x_d = pos_x_a - i_dt*0.5*(vel_x_extrapol + vel_x)*inv_earth_radius;
				pos_y_d = pos_y_a - i_dt*0.5*(vel_y_extrapol + vel_y)*inv_earth_radius;
				pos_z_d = pos_z_a - i_dt*0.5*(vel_z_extrapol + vel_z)*inv_earth_radius;

				ScalarDataArray norm = (pos_x_d*pos_x_d + pos_y_d*pos_y_d + pos_z_d*pos_z_d).inv_sqrt();

				ScalarDataArray new_pos_x_d = pos_x_d*norm;
				ScalarDataArray new_pos_y_d = pos_y_d*norm;
				ScalarDataArray new_pos_z_d = pos_z_d*norm;

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
#endif



#if 0
	/**
	 * Stable extrapolation Two-Time-Level Scheme, Mariano Hortal,
	 *     Development and testing of a new two-time-level semi-lagrangian scheme (settls) in the ECMWF forecast model.
	 * Quaterly Journal of the Royal Meterological Society
	 *
	 * r_d = r_a - dt/2 * (2 * v_n(r_d) - v_{n-1}(r_d) + v_n(r_a))
	 *
	 * v^{iter} := (dt*v_n - dt*0.5*v_{n-1})
	 * r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
	 */
	void semi_lag_departure_points_settls(
			const SphereDataPhysical &i_u_prev,	// Velocities at time t-1
			const SphereDataPhysical &i_v_prev,

			const SphereDataPhysical &i_u, 		// Velocities at time t
			const SphereDataPhysical &i_v,

			const ScalarDataArray &i_posx_a,	// Position of arrival points x / y
			const ScalarDataArray &i_posy_a,

			double i_dt,				///< time step size
			double i_earth_radius,
			ScalarDataArray &o_posx_d, 	///< Position of departure points x / y
			ScalarDataArray &o_posy_d,
			int i_timestepping_order
	)
	{
		if (i_timestepping_order == 1)
		{
			/*
			 * 1st order accurate implementation
			 */
			sphereCoordMinusSurfaceVector(
							i_posx_a,
							i_posy_a,
#if 1
							sample2D.bilinear_scalar(i_u, i_posx_a, i_posy_a, true)/i_earth_radius,
							sample2D.bilinear_scalar(i_v, i_posx_a, i_posy_a, true)/i_earth_radius,
#else
							u_analytical(i_posx_a, i_posy_a)/i_earth_radius,
							v_analytical(i_posx_a, i_posy_a)/i_earth_radius,
#endif
							i_dt,
							o_posx_d,
							o_posy_d
						);
			return;
		}

		/*
		 * TODO:
		 * TODO: Replace this with cartesian coordinates!!!
		 * TODO:
		 */
		std::size_t num_points = i_posx_a.number_of_elements;

		ScalarDataArray u = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_u);
		ScalarDataArray v = Convert_SphereDataPhysical_To_ScalarDataArray::physical_convert(i_v);

		// local dt
		double dt = i_dt;

		// Extrapolate velocities at departure points
		SphereData u_extrapol = 2.0*i_u - i_u_prev;
		SphereData v_extrapol = 2.0*i_v - i_v_prev;

		// Departure point tmp
		ScalarDataArray rx_d_new(num_points);
		ScalarDataArray ry_d_new(num_points);

		// Previous departure point
		ScalarDataArray rx_d_prev = i_posx_a;
		ScalarDataArray ry_d_prev = i_posy_a;

		// initialize departure points with arrival points
		o_posx_d = i_posx_a;
		o_posy_d = i_posy_a;

		double vel_scaling = 1.0/i_earth_radius;

		double diff = -1;
		int iters = 0;
		for (; iters < 100; iters++)
		{
			sphereCoordMinusSurfaceVector(
					i_posx_a,
					i_posy_a,

#if 1
					(0.5 * (u + sample2D.bilinear_scalar(u_extrapol, o_posx_d, o_posy_d, true)))*vel_scaling,
					(0.5 * (v + sample2D.bilinear_scalar(v_extrapol, o_posx_d, o_posy_d, true)))*vel_scaling,
#elif 1
					(0.5 * (u + u_analytical(o_posx_d, o_posy_d)))*vel_scaling,
					(0.5 * (v + v_analytical(o_posx_d, o_posy_d)))*vel_scaling,
#elif 1
					u_analytical(o_posx_d, o_posy_d)*vel_scaling,
					v_analytical(o_posx_d, o_posy_d)*vel_scaling,
#endif

					dt,
					rx_d_new,
					ry_d_new
				);

			/*
			 * TODO: This only works if this is not a point close to the poles
			 */
			diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
			for (std::size_t i = 0; i < num_points; i++)
			{
				o_posx_d.scalar_data[i] = rx_d_new.scalar_data[i];
				o_posy_d.scalar_data[i] = ry_d_new.scalar_data[i];

				// posx \in [-pi/2;pi/2]
				if (o_posy_d.scalar_data[i] > M_PI*0.5)
				{
					o_posy_d.scalar_data[i] = M_PI - o_posy_d.scalar_data[i];
					o_posx_d.scalar_data[i] += M_PI;
				}
				else if (o_posy_d.scalar_data[i] < -M_PI*0.5)
				{
					o_posy_d.scalar_data[i] = -M_PI - o_posy_d.scalar_data[i];
					o_posx_d.scalar_data[i] += M_PI;
				}

				// posx \in [0;2*pi]
				o_posx_d.scalar_data[i] = SphereDataSampler::wrapPeriodic(o_posx_d.scalar_data[i], 2.0*M_PI);

				assert(o_posx_d.scalar_data[i] >= 0);
				assert(o_posx_d.scalar_data[i] < M_PI*2.0);

#if SWEET_DEBUG
				if (o_posy_d.scalar_data[i] < -M_PI*0.5)
					std::cout << o_posy_d.scalar_data[i] << " >= -M_PI*0.5" << std::endl;
				if (o_posy_d.scalar_data[i] > M_PI*0.5)
					std::cout << o_posy_d.scalar_data[i] << " <= M_PI*0.5" << std::endl;
#endif

				assert(o_posy_d.scalar_data[i] >= -M_PI*0.5);
				assert(o_posy_d.scalar_data[i] <= M_PI*0.5);
			}

			if (diff < 1e-10)
			   break;

			rx_d_prev = o_posx_d;
			ry_d_prev = o_posy_d;
		}


		if (iters >= 5)
		{
			std::cout << "WARNING: Too many iterations for SL scheme" << std::endl;
			std::cout << "Iters: " << iters << std::endl;
			std::cout << "DIFF: " << diff << std::endl;
		}

	}
#endif

};

#endif /* SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_ */
