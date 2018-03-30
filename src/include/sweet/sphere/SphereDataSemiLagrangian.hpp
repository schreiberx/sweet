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
	SphereDataConfig *sphereDataConfig;

public:
	SphereDataSemiLagrangian()	:
		sphereDataConfig(nullptr)
	{
	}


	void setup(
		double i_domain_size[2],
		SphereDataConfig *i_sphereDataConfig
	)
	{
		sphereDataConfig = i_sphereDataConfig;
		sample2D.setup(i_domain_size, sphereDataConfig);
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
	void compute_departure_points_settls(
			SphereDataPhysical* i_velocity_field_t_prev[2],	///< velocity field at time n-1
			SphereDataPhysical* i_velocity_field_t[2],		///< velocity field at time n

			ScalarDataArray* i_pos_arrival[2],		///< position at time n+1
			double i_dt,							///< time step size

			ScalarDataArray* o_pos_departure[2],	///< departure points at time n,
			double *i_staggering = nullptr			///< staggering (ux, uy, vx, vy)
	)
	{
		if (i_staggering == nullptr)
		{
			static double constzerostuff[4] = {0,0,0,0};
			i_staggering = constzerostuff;
		}


		std::size_t num_points = i_pos_arrival[0]->number_of_elements;

		/**
		 * Convert velocity components at previous departure points
		 * to ScalarDataArray
		 */

		ScalarDataArray vx_n_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(*i_velocity_field_t_prev[0]);
		ScalarDataArray vy_n_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(*i_velocity_field_t_prev[1]);

		ScalarDataArray vx_n = Convert_SphereData_To_ScalarDataArray::physical_convert(*i_velocity_field_t[0]);
		ScalarDataArray vy_n = Convert_SphereData_To_ScalarDataArray::physical_convert(*i_velocity_field_t[1]);

		ScalarDataArray &rx_a = *i_pos_arrival[0];
		ScalarDataArray &ry_a = *i_pos_arrival[1];

		ScalarDataArray &rx_d = *o_pos_departure[0];
		ScalarDataArray &ry_d = *o_pos_departure[1];


		double dt = i_dt;

		ScalarDataArray vx_iter(num_points);
		ScalarDataArray vy_iter(num_points);

		vx_iter = dt * vx_n - dt*0.5 * vx_n_prev;
		vy_iter = dt * vy_n - dt*0.5 * vy_n_prev;

		ScalarDataArray rx_d_new(num_points);
		ScalarDataArray ry_d_new(num_points);

		ScalarDataArray rx_d_prev = rx_a;
		ScalarDataArray ry_d_prev = ry_a;

		//SphereData* r_d[2] = {&rx_d, &ry_d};

		// initialize departure points with arrival points
		rx_d = rx_a;
		ry_d = ry_a;

		int iters;
		for (iters = 0; iters < 10; iters++)
		{
			// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
			rx_d_new = rx_a - dt*0.5 * vx_n - sample2D.bilinear_scalar(
					Convert_ScalarDataArray_to_SphereData::convert(vx_iter, sphereDataConfig),
					rx_d, ry_d,
					i_staggering[0], i_staggering[1]
			);

			ry_d_new = ry_a - dt*0.5 * vy_n - sample2D.bilinear_scalar(
					Convert_ScalarDataArray_to_SphereData::convert(vy_iter, sphereDataConfig),
					rx_d, ry_d, i_staggering[2], i_staggering[3]);

			double diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();
			rx_d_prev = rx_d_new;
			ry_d_prev = ry_d_new;

			for (	std::size_t i = 0;
					i < num_points;
					i++
			)
			{
				rx_d.scalar_data[i] = sample2D.wrapPeriodic(rx_d_new.scalar_data[i], sample2D.domain_size[0]);
				ry_d.scalar_data[i] = sample2D.wrapPeriodic(ry_d_new.scalar_data[i], sample2D.domain_size[1]);
			}

			if (diff < 1e-8)
				break;
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
			ScalarDataArray &o_posx_d, 	///< Position of departure points x / y
			ScalarDataArray &o_posy_d,

			const SphereStaggering &i_staggering	///< staggering, if any (ux, uy, vx, vy)
	)
	{
		std::size_t num_points = i_posx_a.number_of_elements;

		ScalarDataArray u_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(i_u_prev, false);
		ScalarDataArray v_prev = Convert_SphereData_To_ScalarDataArray::physical_convert(i_v_prev, false);

		ScalarDataArray u = Convert_SphereData_To_ScalarDataArray::physical_convert(i_u, false);
		ScalarDataArray v = Convert_SphereData_To_ScalarDataArray::physical_convert(i_v, false);

		//local dt
		double dt = i_dt;

		//Velocity for iterations
		ScalarDataArray u_iter = dt * u - dt*0.5 * u_prev;
		ScalarDataArray v_iter = dt * v - dt*0.5 * v_prev;

		//Departure point tmp
		ScalarDataArray rx_d_new(num_points);
		ScalarDataArray ry_d_new(num_points);

		//Previous departure point
		ScalarDataArray rx_d_prev = i_posx_a;
		ScalarDataArray ry_d_prev = i_posy_a;

		// initialize departure points with arrival points
		o_posx_d = i_posx_a;
		o_posy_d = i_posy_a;

		int iters = 0;
		for (; iters < 10; iters++)
		{
			//std::cout<<iters<<std::endl;
			// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
			rx_d_new = i_posx_a - dt*0.5 * u - sample2D.bilinear_scalar(
					Convert_ScalarDataArray_to_SphereData::convert(u_iter, sphereDataConfig),
					o_posx_d, o_posy_d, i_staggering.u[0], i_staggering.u[1]
			);
			ry_d_new = i_posy_a - dt*0.5 * v - sample2D.bilinear_scalar(
					Convert_ScalarDataArray_to_SphereData::convert(v_iter, sphereDataConfig),
					o_posx_d, o_posy_d, i_staggering.v[0], i_staggering.v[1]
			);

			double diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();
			rx_d_prev = rx_d_new;
			ry_d_prev = ry_d_new;

			for (std::size_t i = 0; i < num_points; i++)
			{
				o_posx_d.scalar_data[i] = sample2D.wrapPeriodic(rx_d_new.scalar_data[i], sample2D.domain_size[0]);
				o_posy_d.scalar_data[i] = sample2D.wrapPeriodic(ry_d_new.scalar_data[i], sample2D.domain_size[1]);
			}

			if (diff < 1e-8)
			   break;
		}
	}
};

#endif /* SRC_INCLUDE_SWEET_SPHEREDATASEMILAGRANGIAN_HPP_ */
