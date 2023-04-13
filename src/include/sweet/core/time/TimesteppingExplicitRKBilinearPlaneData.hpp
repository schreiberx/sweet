
#ifndef PLANEDATA_TIMESTEPPING_EXPLICIT_RK_HPP__
#define PLANEDATA_TIMESTEPPING_EXPLICIT_RK_HPP__

#include "../plane/PlaneData_Spectral.hpp"

namespace sweet
{

class TimesteppingExplicitRKBilinearPlaneData
{
	// runge kutta data storages
	PlaneData_Spectral** RK_h_t;
	PlaneData_Spectral** RK_u_t;
	PlaneData_Spectral** RK_v_t;

	int runge_kutta_order;

public:
	TimesteppingExplicitRKBilinearPlaneData()	:
		RK_h_t(nullptr),
		RK_u_t(nullptr),
		RK_v_t(nullptr),
		runge_kutta_order(-1)
	{
	}



	void setupBuffers(
			const PlaneData_Config *i_planeDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
	)
	{
		if (RK_h_t != nullptr)	///< already allocated?
			return;

		runge_kutta_order = i_rk_order;
		int N = i_rk_order;

		if (N <= 0 || N > 4)
			SWEETError("Invalid order for RK time stepping (Please set --timestepping-order and/or --timestepping-order2)");

		RK_h_t = new PlaneData_Spectral*[N];
		RK_u_t = new PlaneData_Spectral*[N];
		RK_v_t = new PlaneData_Spectral*[N];

		for (int i = 0; i < N; i++)
		{
			RK_h_t[i] = new PlaneData_Spectral(i_planeDataConfig);
			RK_u_t[i] = new PlaneData_Spectral(i_planeDataConfig);
			RK_v_t[i] = new PlaneData_Spectral(i_planeDataConfig);
		}
	}



	~TimesteppingExplicitRKBilinearPlaneData()
	{
		int N = runge_kutta_order;

		if (RK_h_t != nullptr)
		{
			for (int i = 0; i < N; i++)
			{
				delete RK_h_t[i];
				delete RK_u_t[i];
				delete RK_v_t[i];
			}

			delete [] RK_h_t;
			delete [] RK_u_t;
			delete [] RK_v_t;

			RK_h_t = nullptr;
			RK_u_t = nullptr;
			RK_v_t = nullptr;
		}
	}



	/**
	 * execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void runTimestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const PlaneData_Spectral &i_P_A,	///< prognostic variables
					const PlaneData_Spectral &i_u_A,	///< prognostic variables
					const PlaneData_Spectral &i_v_A,	///< prognostic variables

					const PlaneData_Spectral &i_P_B,	///< prognostic variables
					const PlaneData_Spectral &i_u_B,	///< prognostic variables
					const PlaneData_Spectral &i_v_B,	///< prognostic variables

					PlaneData_Spectral &o_P_t,	///< time updates
					PlaneData_Spectral &o_u_t,	///< time updates
					PlaneData_Spectral &o_v_t,	///< time updates

					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			PlaneData_Spectral &i_var0_A,
			PlaneData_Spectral &i_var1_A,
			PlaneData_Spectral &i_var2_A,

			PlaneData_Spectral &i_var0_B,
			PlaneData_Spectral &i_var1_B,
			PlaneData_Spectral &i_var2_B,

			PlaneData_Spectral &o_var0,
			PlaneData_Spectral &o_var1,
			PlaneData_Spectral &o_var2,

			double i_dt = 0,				///< Use this time step size
			int i_runge_kutta_order = 1,	///< Order of RK time stepping
			double i_simulation_time = -1	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		setupBuffers(i_var0_A.planeDataConfig, i_runge_kutta_order);

		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A,	// input
					i_var1_A,
					i_var2_A,
					i_var0_B,	// input
					i_var1_B,
					i_var2_B,
					*RK_h_t[0],	// output
					*RK_u_t[0],
					*RK_v_t[0],
					i_simulation_time
			);

			o_var0 = i_dt**RK_h_t[0];
			o_var1 = i_dt**RK_u_t[0];
			o_var2 = i_dt**RK_v_t[0];
		}
		else if (i_runge_kutta_order == 2)
		{
			// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
			// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
			/*
			 * c     a
			 * 0   |
			 * 1/2 | 1/2
			 * --------------
			 *     | 0   1    b
			 */
			double a2[1] = {0.5};
			double b[2] = {0.0, 1.0};
			double c[1] = {0.5};

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A,
					i_var1_A,
					i_var2_A,
					i_var0_B,
					i_var1_B,
					i_var2_B,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + ( i_dt*a2[0]*(*RK_h_t[0]) ),
					i_var1_A + ( i_dt*a2[0]*(*RK_u_t[0]) ),
					i_var2_A + ( i_dt*a2[0]*(*RK_v_t[0]) ),
					i_var0_B + ( i_dt*a2[0]*(*RK_h_t[0]) ),
					i_var1_B + ( i_dt*a2[0]*(*RK_u_t[0]) ),
					i_var2_B + ( i_dt*a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					i_simulation_time + c[0]*i_dt
			);

			o_var0 = i_dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_h_t[1]) );
			o_var1 = i_dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_u_t[1]) );
			o_var2 = i_dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_v_t[1]) );
		}
		else if (i_runge_kutta_order == 3)
		{
			// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
			// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
			/*
			 * c     a
			 * 0   |
			 * 1/3 | 1/3
			 * 2/3 | 0    2/3
			 * --------------
			 *     | 1/4  0   3/4
			 */
			double a2[1] = {1.0/3.0};
			double a3[2] = {0.0, 2.0/3.0};
			double b[3] = {1.0/4.0, 0.0, 3.0/4.0};
			double c[2] = {1.0/3.0, 2.0/3.0};

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A,
					i_var1_A,
					i_var2_A,
					i_var0_B,
					i_var1_B,
					i_var2_B,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + i_dt*( a2[0]*(*RK_h_t[0]) ),
					i_var1_A + i_dt*( a2[0]*(*RK_u_t[0]) ),
					i_var2_A + i_dt*( a2[0]*(*RK_v_t[0]) ),
					i_var0_B + i_dt*( a2[0]*(*RK_h_t[0]) ),
					i_var1_B + i_dt*( a2[0]*(*RK_u_t[0]) ),
					i_var2_B + i_dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + i_dt*( a3[0]*(*RK_h_t[0]) + a3[1]*(*RK_h_t[1]) ),
					i_var1_A + i_dt*( a3[0]*(*RK_u_t[0]) + a3[1]*(*RK_u_t[1]) ),
					i_var2_A + i_dt*( a3[0]*(*RK_v_t[0]) + a3[1]*(*RK_v_t[1]) ),
					i_var0_B + i_dt*( a3[0]*(*RK_h_t[0]) + a3[1]*(*RK_h_t[1]) ),
					i_var1_B + i_dt*( a3[0]*(*RK_u_t[0]) + a3[1]*(*RK_u_t[1]) ),
					i_var2_B + i_dt*( a3[0]*(*RK_v_t[0]) + a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					i_simulation_time + c[1]*i_dt
			);

			o_var0 = i_dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) );
			o_var1 = i_dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) );
			o_var2 = i_dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) );
		}
		else if (i_runge_kutta_order == 4)
		{
			// See https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods
			// See https://de.wikipedia.org/wiki/Runge-Kutta-Verfahren
			/*
			 * c     a
			 * 0   |
			 * 1/2 | 1/2
			 * 1/2 | 0    1/2
			 * 1   | 0    0    1
			 * --------------
			 *     | 1/6  1/3  1/3  1/6
			 */
			double a2[1] = {0.5};
			double a3[2] = {0.0, 0.5};
			double a4[3] = {0.0, 0.0, 1.0};
			double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
			double c[3] = {0.5, 0.5, 1.0};

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A,
					i_var1_A,
					i_var2_A,
					i_var0_B,
					i_var1_B,
					i_var2_B,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + i_dt*( a2[0]*(*RK_h_t[0]) ),
					i_var1_A + i_dt*( a2[0]*(*RK_u_t[0]) ),
					i_var2_A + i_dt*( a2[0]*(*RK_v_t[0]) ),
					i_var0_B + i_dt*( a2[0]*(*RK_h_t[0]) ),
					i_var1_B + i_dt*( a2[0]*(*RK_u_t[0]) ),
					i_var2_B + i_dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + i_dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_h_t[1]) ),
					i_var1_A + i_dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_u_t[1]) ),
					i_var2_A + i_dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_v_t[1]) ),
					i_var0_B + i_dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_h_t[1]) ),
					i_var1_B + i_dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_u_t[1]) ),
					i_var2_B + i_dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					i_simulation_time + c[1]*i_dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					i_var0_A + i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_t[2]) ),
					i_var1_A + i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_t[2]) ),
					i_var2_A + i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_t[2]) ),
					i_var0_B + i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_t[2]) ),
					i_var1_B + i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_t[2]) ),
					i_var2_B + i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_t[2]) ),
					*RK_h_t[3],
					*RK_u_t[3],
					*RK_v_t[3],
					i_simulation_time + c[2]*i_dt
			);

			o_var0 = i_dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) + (b[3]*(*RK_h_t[3])) );
			o_var1 = i_dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) + (b[3]*(*RK_u_t[3])) );
			o_var2 = i_dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) + (b[3]*(*RK_v_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}
};

}

#endif
