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
		sample2D.setup(i_domain_size, sphereDataConfig);
	}


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



	inline
	static
	void cartToAngleCoord(
			const double i_x[3],
			double *o_lon,
			double *o_lat
	)
	{
		*o_lon = std::acos(i_x[0]/std::sqrt(i_x[0]*i_x[0] + i_x[1]*i_x[1]));
		if (i_x[1] < 0)
			*o_lon = 2.0*M_PI-*o_lon;

		*o_lat = std::acos(-i_x[2]);
		*o_lat -= M_PI*0.5;
	}


	inline
	static
	double length(
			const double *i_x
	)
	{
		return std::sqrt(i_x[0]*i_x[0] + i_x[1]*i_x[1] + i_x[2]*i_x[2]);
	}


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


	/**
	 * Add a vector on the sphere to a position
	 */
	void sphereCoordPlusSurfaceVector(
			const ScalarDataArray &i_pos_lon,	/// longitude angle
			const ScalarDataArray &i_pos_lat,	/// latitude angle

			const ScalarDataArray &i_vec_lon,	/// velocity along longitude
			const ScalarDataArray &i_vec_lat,	/// velocity along latitude

			ScalarDataArray &o_pos_lon,
			ScalarDataArray &o_pos_lat
		)
	{
		// TODO: make this vectorizable

		// get lat-lon tangential basis
#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
		for (std::size_t i = 0; i < i_pos_lon.number_of_elements; i++)
		{
			// compute 3D pos
			double tanBasis[9];
			angleToTangentialSpace(i_pos_lon.scalar_data[i], i_pos_lat.scalar_data[i], tanBasis);

			// compute 3D velocity
			double velVec[3];
			for (int j = 0; j < 3; j++)
				velVec[j] = tanBasis[j+3]*i_vec_lon.scalar_data[i] + tanBasis[j+6]*i_vec_lat.scalar_data[i];

#if 1
			// velocity basis
			double velBasis[9];
			for (int j = 0; j < 3; j++)
				velBasis[j] = tanBasis[j];

			double vel_length = std::sqrt(i_vec_lon.scalar_data[i]*i_vec_lon.scalar_data[i] + i_vec_lat.scalar_data[i]*i_vec_lat.scalar_data[i]);
			double inv_vel_len = 1.0/vel_length;
			for (int j = 0; j < 3; j++)
				velBasis[3+j] = velVec[j]*inv_vel_len;

			cross(&velBasis[0], &velBasis[3], &velBasis[6]);

/*
			for (int j = 0; j < 3; j++)
			{
				std::cout << velBasis[3*j+0] << ", " << velBasis[3*j+1] << ", " << velBasis[3*j+2] << std::endl;
			}
			std::cout << std::endl;
*/

			// Follow vel along unit circle
			// start at relative position (1,0)
			double pos[2];
			pos[0] = std::cos(vel_length)-1.0;
			pos[1] = std::sin(vel_length);

			// update velVec
#if 1
			for (int j = 0; j < 3; j++)
			{
				velVec[j] = 0;
				for (int i = 0; i < 2; i++)
				{
					// 2D velocity basis is made up by 2nd and 3rd column
					velVec[j] += velBasis[j+i*3]*pos[i];
				}
			}
#else
			velVec[0] = 0;
			for (int i = 0; i < 2; i++)
				velVec[0] += velBasis[i*3]*pos[i];

			velVec[1] = 0;
			for (int i = 0; i < 2; i++)
				velVec[1] += velBasis[3+i*3]*pos[i];
#endif

			// add
			double finalPos[3];
			for (int j = 0; j < 3; j++)
				finalPos[j] = tanBasis[j] + velVec[j];

//			std::cout << tanBasis[0] << ", " << tanBasis[1] << ", " << tanBasis[2] << std::endl;
//			std::cout << velVec[0] << ", " << velVec[1] << ", " << velVec[2] << std::endl;
//			std::cout << length(finalPos) << std::endl;

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
		}
	}


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
			ScalarDataArray &o_posy_d
	)
	{
		std::size_t num_points = i_posx_a.number_of_elements;

		ScalarDataArray u_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(i_u_prev, false);
		ScalarDataArray v_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(i_v_prev, false);

		ScalarDataArray u = Convert_SphereData_To_ScalarDataArray::physical_convert(i_u, false);
		ScalarDataArray v = Convert_SphereData_To_ScalarDataArray::physical_convert(i_v, false);

		// local dt
		double dt = i_dt;

		// Velocity for iterations
		SphereData u_iter = Convert_ScalarDataArray_to_SphereData::convert(dt * u - dt*0.5 * u_prev, sphereDataConfig);
		SphereData v_iter = Convert_ScalarDataArray_to_SphereData::convert(dt * v - dt*0.5 * v_prev, sphereDataConfig);

		u_iter.physical_set_zero();
		v_iter.physical_set_zero();


		// Departure point tmp
		ScalarDataArray rx_d_new(num_points);
		ScalarDataArray ry_d_new(num_points);

		// Previous departure point
		ScalarDataArray rx_d_prev = i_posx_a;
		ScalarDataArray ry_d_prev = i_posy_a;

		// initialize departure points with arrival points
		o_posx_d = i_posx_a;
		o_posy_d = i_posy_a;

		std::cout << "TODO: Figure out where the scaling 2 factor comes from!" << std::endl;
		double vel_scaling = 2.0/i_earth_radius;

		int iters = 0;
		for (; iters < 10; iters++)
		{
			sphereCoordPlusSurfaceVector(
					i_posx_a,
					i_posy_a,
					(-dt*0.5 * u - sample2D.bilinear_scalar(u_iter, o_posx_d, o_posy_d))*vel_scaling,
					(-dt*0.5 * v - sample2D.bilinear_scalar(v_iter, o_posx_d, o_posy_d))*vel_scaling,
					rx_d_new,
					ry_d_new
				);

			double diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();

#if SWEET_SPACE_THREADING
#pragma omp parallel for
#endif
			for (std::size_t i = 0; i < num_points; i++)
			{
				// posx \in [0;2*pi]
				o_posx_d.scalar_data[i] = SphereDataSampler::wrapPeriodic(rx_d_new.scalar_data[i], 2.0*M_PI);
				assert(o_posx_d.scalar_data[i] >= 0);
				assert(o_posx_d.scalar_data[i] < M_PI*2.0);

				// posx \in [-pi/2;pi/2]
				o_posy_d.scalar_data[i] = ry_d_new.scalar_data[i];
				if (o_posy_d.scalar_data[i] > M_PI*0.5)
					o_posy_d.scalar_data[i] = M_PI - o_posy_d.scalar_data[i];
				else if (o_posy_d.scalar_data[i] < -M_PI*0.5)
					o_posy_d.scalar_data[i] = -M_PI - o_posy_d.scalar_data[i];

				assert(o_posy_d.scalar_data[i] >= -M_PI*0.5);
				assert(o_posy_d.scalar_data[i] <= M_PI*0.5);
			}

			if (diff < 1e-8)
			   break;

			rx_d_prev = o_posx_d;
			ry_d_prev = o_posy_d;
		}


		if (iters == 10)
			std::cout << "WARNING: Too many iterations for SL scheme" << std::endl;
	}
};

#endif /* SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_ */
