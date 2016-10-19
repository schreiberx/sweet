/*
 * SemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_
#define SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_

#include "PlaneData.hpp"
#include "PlaneDataSampler.hpp"

class SemiLagrangian
{
	PlaneDataSampler sample2D;
	PlaneDataConfig *planeDataConfig;

public:
	SemiLagrangian()	:
		planeDataConfig(nullptr)
	{
	}


	void setup(
		double i_domain_size[2],
		PlaneDataConfig *i_planeDataConfig
	)
	{
		planeDataConfig = i_planeDataConfig;
		sample2D.setup(i_domain_size, planeDataConfig);
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
			PlaneData* i_velocity_field_t_prev[2],	///< velocity field at time n-1
			PlaneData* i_velocity_field_t[2],	///< velocity field at time n
			PlaneData* i_pos_arrival[2],			///< position at time n+1
			double i_dt,							///< time step size
			PlaneData* o_pos_departure[2],			///< departure points at time n,
			double *i_staggering = nullptr				///< staggering (ux, uy, vx, vy)
	)
	{
		if (i_staggering == nullptr)
		{
			static double constzerostuff[4] = {0,0,0,0};
			i_staggering = constzerostuff;
		}

		PlaneData &vx_n_prev = *i_velocity_field_t_prev[0];
		PlaneData &vy_n_prev = *i_velocity_field_t_prev[1];

		PlaneData &vx_n = *i_velocity_field_t[0];
		PlaneData &vy_n = *i_velocity_field_t[1];

		PlaneData &rx_a = *i_pos_arrival[0];
		PlaneData &ry_a = *i_pos_arrival[1];

		PlaneData &rx_d = *o_pos_departure[0];
		PlaneData &ry_d = *o_pos_departure[1];


		double dt = i_dt;

		PlaneData vx_iter(vx_n_prev.planeDataConfig);
		PlaneData vy_iter(vx_n_prev.planeDataConfig);

		vx_iter = dt * vx_n - dt*0.5 * vx_n_prev;
		vy_iter = dt * vy_n - dt*0.5 * vy_n_prev;

		PlaneData rx_d_new(vx_n_prev.planeDataConfig);
		PlaneData ry_d_new(vx_n_prev.planeDataConfig);

		PlaneData rx_d_prev = rx_a;
		PlaneData ry_d_prev = ry_a;	// TODO: is rx_a correct or should it be ry_a?

		//PlaneData* r_d[2] = {&rx_d, &ry_d};

		// initialize departure points with arrival points
		rx_d = rx_a;
		ry_d = ry_a;

		int iters;
		for (iters = 0; iters < 10; iters++)
		{
			// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
			rx_d_new = rx_a - dt*0.5 * vx_n - sample2D.bilinear_scalar(vx_iter, rx_d, ry_d, i_staggering[0], i_staggering[1]);
			ry_d_new = ry_a - dt*0.5 * vy_n - sample2D.bilinear_scalar(vy_iter, rx_d, ry_d, i_staggering[2], i_staggering[3]);

			//std::cout << "WHATS GOING ON HERE?!?" << std::endl;
//			std::cout << rx_d_new << std::endl;
//			exit(1);

			double diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();
			rx_d_prev = rx_d_new;
			ry_d_prev = ry_d_new;

			for (	std::size_t i = 0;
					i < rx_d.planeDataConfig->physical_data_size[0]*rx_d.planeDataConfig->physical_data_size[1];
					i++
			)
			{
				rx_d.physical_space_data[i] = sample2D.wrapPeriodic(rx_d_new.physical_space_data[i], sample2D.domain_size[0]);
				ry_d.physical_space_data[i] = sample2D.wrapPeriodic(ry_d_new.physical_space_data[i], sample2D.domain_size[1]);
			}

			if (diff < 1e-8)
				break;
//			std::cout << iters << ": " << diff << std::endl;
		}
//		std::cout << iters << std::endl;
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
			PlaneData &i_u_prev, // Velocities at time t-1
			PlaneData &i_v_prev,
			PlaneData &i_u, // Velocities at time t
			PlaneData &i_v,
			PlaneData &i_posx_a, // Position of arrival points x / y
			PlaneData &i_posy_a,
			double i_dt,			///< time step size
			PlaneData &o_posx_d, // Position of departure points x / y
			PlaneData &o_posy_d,
			double *i_staggering = nullptr	///< staggering, if any (ux, uy, vx, vy)
	)
	{
		i_u_prev.request_data_physical();
		i_v_prev.request_data_physical();
		i_u.request_data_physical();
		i_v.request_data_physical();
		i_posx_a.request_data_physical();
		i_posy_a.request_data_physical();

		if (i_staggering == nullptr)
		{
			static double constzerostuff[4] = {0,0,0,0};
			i_staggering = constzerostuff;
		}

		//local dt
		double dt = i_dt;

		//Velocity for iterations
		PlaneData u_iter(i_u.planeDataConfig);
		PlaneData v_iter(i_v.planeDataConfig);

		//Time Extrapolation
		u_iter = dt * i_u - dt*0.5 * i_u_prev;
		v_iter = dt * i_v - dt*0.5 * i_v_prev;

		//Departure point tmp
		PlaneData rx_d_new(i_u.planeDataConfig);
		PlaneData ry_d_new(i_v.planeDataConfig);

		//Previous departure point
		PlaneData rx_d_prev = i_posx_a;
		PlaneData ry_d_prev = i_posy_a;

		// initialize departure points with arrival points
		o_posx_d = i_posx_a;
		o_posy_d = i_posy_a;

		//std::cout<<"u"<<std::endl;
		//i_u.printArrayData();
		//std::cout<<"v"<<std::endl;
		//i_v.printArrayData();

		//std::cout<<"i_posx_a"<<std::endl;
		//i_posx_a.printArrayData();
		//std::cout<<"i_posy_a"<<std::endl;
		//i_posy_a.printArrayData();

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		o_posx_d.physical_space_data_valid = true;
		o_posy_d.physical_space_data_valid = true;
		o_posx_d.spectral_space_data_valid = false;
		o_posy_d.spectral_space_data_valid = false;

		assert(i_posx_a.physical_space_data_valid);
		assert(i_posy_a.physical_space_data_valid);

		assert(i_u.physical_space_data_valid);
		assert(i_v.physical_space_data_valid);

		assert(i_u_prev.physical_space_data_valid);
		assert(i_v_prev.physical_space_data_valid);
#endif

		int iters = 0;
		for (; iters < 10; iters++)
		{
			//std::cout<<iters<<std::endl;
			// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
			rx_d_new = i_posx_a - dt*0.5 * i_u - sample2D.bilinear_scalar(u_iter, o_posx_d, o_posy_d, i_staggering[0], i_staggering[1]);
			ry_d_new = i_posy_a - dt*0.5 * i_v - sample2D.bilinear_scalar(v_iter, o_posx_d, o_posy_d, i_staggering[2], i_staggering[3]);

			//std::cout<<"u_iter"<<std::endl;
			//u_iter.printArrayData();
			//std::cout<<"v_iter"<<std::endl;
			//v_iter.printArrayData();

			double diff = (rx_d_new - rx_d_prev).reduce_maxAbs() + (ry_d_new - ry_d_prev).reduce_maxAbs();
			rx_d_prev = rx_d_new;
			ry_d_prev = ry_d_new;

			rx_d_new.request_data_physical();
			ry_d_new.request_data_physical();


			for (std::size_t i = 0; i < o_posx_d.planeDataConfig->physical_array_data_number_of_elements; i++)
			{
				o_posx_d.physical_space_data[i] = sample2D.wrapPeriodic(rx_d_new.physical_space_data[i], sample2D.domain_size[0]);
				 //std::cout<<i<<'\t'<<o_posx_d.array_data_cartesian_space[i]<<'\t'<< rx_d_new.array_data_cartesian_space[i] << '\t'<< sample2D.domain_size[0] <<std::endl;
				o_posy_d.physical_space_data[i] = sample2D.wrapPeriodic(ry_d_new.physical_space_data[i], sample2D.domain_size[1]);
			}

			//std::cout<<"Departure_x"<<std::endl;
			//o_posx_d.printArrayData();
			//std::cout<<"Departure_y"<<std::endl;
			//o_posy_d.printArrayData();


			if (diff < 1e-8)
			   break;

			//std::cout << iters << " : " << diff << std::endl;

		}
		//std::cout << iters << std::endl;
	}
};

#endif /* SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_ */
