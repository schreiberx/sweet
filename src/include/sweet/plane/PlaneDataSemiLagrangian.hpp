/*
 * SemiLangrangian.hpp
 *
 *  Created on: 5 Dec 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_
#define SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_

#include <sweet/plane/Convert_PlaneData_to_ScalarDataArray.hpp>
#include <sweet/plane/Convert_ScalarDataArray_to_PlaneData.hpp>
#include <sweet/plane/PlaneStaggering.hpp>
#include "PlaneData.hpp"
#include "PlaneDataSampler.hpp"
#include <sweet/ScalarDataArray.hpp>
#include <sweet/SimulationBenchmarkTiming.hpp>

class PlaneDataSemiLagrangian
{
	PlaneDataSampler sample2D;
	const PlaneDataConfig *planeDataConfig;

public:
	PlaneDataSemiLagrangian()	:
		planeDataConfig(nullptr)
	{
	}


	void setup(
		double i_domain_size[2],
		const PlaneDataConfig *i_planeDataConfig
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
	void semi_lag_departure_points_settls(
			const PlaneData &i_u_prev,	// Velocities at time t-1
			const PlaneData &i_v_prev,
			const PlaneData &i_u, 		// Velocities at time t
			const PlaneData &i_v,

			const ScalarDataArray &i_posx_a,	// Position of arrival points x / y
			const ScalarDataArray &i_posy_a,

			double i_dt,				///< time step size
			ScalarDataArray &o_posx_d, 	///< Position of departure points x / y
			ScalarDataArray &o_posy_d,

			double i_domain_size[2],	///< domain size

			const Staggering *i_staggering,	///< staggering, if any (ux, uy, vx, vy)

			int i_timestepping_order,

			int max_iters,
			double convergence_tolerance
	)
	{
#if SWEET_BENCHMARK_TIMINGS
	SimulationBenchmarkTimings::getInstance().main_timestepping_semi_lagrangian.start();
#endif

		Staggering s;
		if (i_staggering == nullptr)
		{
			s.setup_a_staggering();
			i_staggering = &(const Staggering&)s;
		}

		std::size_t num_points = i_posx_a.number_of_elements;

		if (i_timestepping_order == 1)
		{
			o_posx_d = i_posx_a - i_dt*Convert_PlaneData_To_ScalarDataArray::physical_convert(i_u, false);
			o_posy_d = i_posy_a - i_dt*Convert_PlaneData_To_ScalarDataArray::physical_convert(i_v, false);
		}
		else if (i_timestepping_order == 2)
		{
			ScalarDataArray u_prev = Convert_PlaneData_To_ScalarDataArray::physical_convert(i_u_prev, false);
			ScalarDataArray v_prev = Convert_PlaneData_To_ScalarDataArray::physical_convert(i_v_prev, false);

			ScalarDataArray u = Convert_PlaneData_To_ScalarDataArray::physical_convert(i_u, false);
			ScalarDataArray v = Convert_PlaneData_To_ScalarDataArray::physical_convert(i_v, false);

			double dt = i_dt;

			//Velocity for iterations
			//ScalarDataArray u_iter = dt * u - dt*0.5 * u_prev;
			//ScalarDataArray v_iter = dt * v - dt*0.5 * v_prev;
#if 1
			PlaneData u_extrap = Convert_ScalarDataArray_to_PlaneData::convert(2.0*u - u_prev, planeDataConfig);
			PlaneData v_extrap = Convert_ScalarDataArray_to_PlaneData::convert(2.0*v - v_prev, planeDataConfig);
#else       //To avoid multi-step method
			PlaneData u_extrap = Convert_ScalarDataArray_to_PlaneData::convert(u , planeDataConfig);
			PlaneData v_extrap = Convert_ScalarDataArray_to_PlaneData::convert(v , planeDataConfig);
#endif
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
			double diff = 999;
			for (; iters < max_iters; iters++)
			{
				// r_d = r_a - dt/2 * v_n(r_d) - v^{iter}(r_d)
				rx_d_new = i_posx_a - dt*0.5 *
					(u + sample2D.bilinear_scalar(
						u_extrap,
						o_posx_d, o_posy_d, i_staggering->u[0], i_staggering->u[1]
				));
				ry_d_new = i_posy_a - dt*0.5 *
					(v + sample2D.bilinear_scalar(
						v_extrap,
						o_posx_d, o_posy_d, i_staggering->v[0], i_staggering->v[1]
				));

				diff = (rx_d_new - rx_d_prev).reduce_maxAbs()/i_domain_size[0] + (ry_d_new - ry_d_prev).reduce_maxAbs()/i_domain_size[1];
				rx_d_prev = rx_d_new;
				ry_d_prev = ry_d_new;

				SWEET_THREADING_SPACE_PARALLEL_FOR
				for (std::size_t i = 0; i < num_points; i++)
				{
					o_posx_d.scalar_data[i] = sample2D.wrapPeriodic(rx_d_new.scalar_data[i], sample2D.domain_size[0]);
					o_posy_d.scalar_data[i] = sample2D.wrapPeriodic(ry_d_new.scalar_data[i], sample2D.domain_size[1]);
				}

				if (diff < convergence_tolerance)
				   break;
			}


			if (convergence_tolerance > 0)
			{
				if (diff > convergence_tolerance)
				{
					std::cout << "WARNING: Over convergence tolerance" << std::endl;
					std::cout << "+ Iterations: " << iters << std::endl;
					std::cout << "+ maxAbs: " << diff << std::endl;
					std::cout << "+ Convergence tolerance: " << convergence_tolerance << std::endl;
				}
			}
		}
		else
		{
			FatalError("This time integration order is not implemented");
		}


#if SWEET_BENCHMARK_TIMINGS
		SimulationBenchmarkTimings::getInstance().main_timestepping_semi_lagrangian.stop();
#endif

	}
};

#endif /* SRC_INCLUDE_SWEET_PLANEDATASEMILAGRANGIAN_HPP_ */
