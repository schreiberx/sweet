
#ifndef TIMESTEPPING_RK_HPP
#define TIMESTEPPING_RK_HPP

#include "PlaneData.hpp"

class PlaneDataTimesteppingRK
{
	// runge kutta data storages
	PlaneData** RK_h_t;
	PlaneData** RK_u_t;
	PlaneData** RK_v_t;

	int runge_kutta_order;

public:
	PlaneDataTimesteppingRK()	:
		RK_h_t(nullptr),
		RK_u_t(nullptr),
		RK_v_t(nullptr),
		runge_kutta_order(-1)
	{
	}



	void setupBuffers(
			const PlaneDataConfig *i_planeDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
	)
	{
		if (RK_h_t != nullptr)	///< already allocated?
			return;

		runge_kutta_order = i_rk_order;
		int N = i_rk_order;

		if (N <= 0 || N > 4)
			FatalError("Invalid order for RK time stepping");

		RK_h_t = new PlaneData*[N];
		RK_u_t = new PlaneData*[N];
		RK_v_t = new PlaneData*[N];

		for (int i = 0; i < N; i++)
		{
			RK_h_t[i] = new PlaneData(i_planeDataConfig);
			RK_u_t[i] = new PlaneData(i_planeDataConfig);
			RK_v_t[i] = new PlaneData(i_planeDataConfig);
		}
	}



	~PlaneDataTimesteppingRK()
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
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const PlaneData &i_P,	///< prognostic variables
					const PlaneData &i_u,	///< prognostic variables
					const PlaneData &i_v,	///< prognostic variables

					PlaneData &o_P_t,	///< time updates
					PlaneData &o_u_t,	///< time updates
					PlaneData &o_v_t,	///< time updates

					double &o_dt,			///< time step restriction
					double i_use_fixed_dt,	///< if this value is not equal to 0,
											///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			PlaneData &io_var0,
			PlaneData &io_var1,
			PlaneData &io_var2,

			double &o_dt,					///< return time step size for the computed time step

			double i_use_fixed_dt = 0,		///< If this value is not equal to 0,
											///< Use this time step size instead of computing one
											///< This also sets o_dt = i_use_fixed_dt

			int i_runge_kutta_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1,	///< Current simulation time.
											///< This gets e.g. important for tidal waves

			double i_max_simulation_time = std::numeric_limits<double>::infinity()	///< limit the maximum simulation time
	)
	{
		setupBuffers(io_var0.planeDataConfig, i_runge_kutta_order);

		double &dt = o_dt;
		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0,	// input
					io_var1,
					io_var2,
					*RK_h_t[0],	// output
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

			io_var0 += dt**RK_h_t[0];
			io_var1 += dt**RK_u_t[0];
			io_var2 += dt**RK_v_t[0];

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

			double dummy_dt = -1;

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0,
					io_var1,
					io_var2,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0 + ( dt*a2[0]*(*RK_h_t[0]) ),
					io_var1 + ( dt*a2[0]*(*RK_u_t[0]) ),
					io_var2 + ( dt*a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			io_var0 += dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_h_t[1]) );
			io_var1 += dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_u_t[1]) );
			io_var2 += dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_v_t[1]) );
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

			double dummy_dt;

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0,
					io_var1,
					io_var2,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0	+ dt*( a2[0]*(*RK_h_t[0]) ),
					io_var1	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_var2	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0	+ dt*( a3[0]*(*RK_h_t[0]) + a3[1]*(*RK_h_t[1]) ),
					io_var1	+ dt*( a3[0]*(*RK_u_t[0]) + a3[1]*(*RK_u_t[1]) ),
					io_var2	+ dt*( a3[0]*(*RK_v_t[0]) + a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			io_var0 += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) );
			io_var1 += dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) );
			io_var2 += dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) );
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

			double dummy_dt;

			// STAGE 1
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0,
					io_var1,
					io_var2,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0	+ dt*( a2[0]*(*RK_h_t[0]) ),
					io_var1	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_var2	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0	+ dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_h_t[1]) ),
					io_var1	+ dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_u_t[1]) ),
					io_var2	+ dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_var0	+ dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_t[2]) ),
					io_var1	+ dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_t[2]) ),
					io_var2	+ dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_t[2]) ),
					*RK_h_t[3],
					*RK_u_t[3],
					*RK_v_t[3],
					dummy_dt,
					dt,
					i_simulation_time + c[2]*dt
			);


			io_var0 += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) + (b[3]*(*RK_h_t[3])) );
			io_var1 += dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) + (b[3]*(*RK_u_t[3])) );
			io_var2 += dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) + (b[3]*(*RK_v_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}
};

#endif
