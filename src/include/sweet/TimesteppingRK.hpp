
#ifndef TIMESTEPPING_RK_HPP
#define TIMESTEPPING_RK_HPP


class TimesteppingRK
{
	// runge kutta data storages
	DataArray<2>** RK_h_t;
	DataArray<2>** RK_u_t;
	DataArray<2>** RK_v_t;

	int runge_kutta_order;

public:
	TimesteppingRK()	:
		RK_h_t(nullptr),
		RK_u_t(nullptr),
		RK_v_t(nullptr),
		runge_kutta_order(-1)
	{
	}



	void setupBuffers(
			const DataArray<2> &i_test_buffer,	///< array of example data to know dimensions of buffers
			int i_rk_order			///< Order of Runge-Kutta method
	)
	{
		if (RK_h_t != nullptr)	///< already allocated?
			return;

		runge_kutta_order = i_rk_order;
		int N = i_rk_order;

		RK_h_t = new DataArray<2>*[N];
		RK_u_t = new DataArray<2>*[N];
		RK_v_t = new DataArray<2>*[N];

		for (int i = 0; i < N; i++)
		{
			RK_h_t[i] = new DataArray<2>(i_test_buffer.resolution);
			RK_u_t[i] = new DataArray<2>(i_test_buffer.resolution);
			RK_v_t[i] = new DataArray<2>(i_test_buffer.resolution);
		}
	}



	~TimesteppingRK()
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



	template <class BaseClass>
	void run_rk_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const DataArray<2> &i_P,	///< prognostic variables
					const DataArray<2> &i_u,	///< prognostic variables
					const DataArray<2> &i_v,	///< prognostic variables

					DataArray<2> &o_P_t,	///< time updates
					DataArray<2> &o_u_t,	///< time updates
					DataArray<2> &o_v_t,	///< time updates

					double &o_dt,			///< time step restriction
					double i_fixed_dt,		///< if this value is not equal to 0,
											///< use this time step size instead of computing one
					double i_simulation_time
			),

			DataArray<2> &io_h,
			DataArray<2> &io_u,
			DataArray<2> &io_v,

			double &o_dt,
			double i_fixed_dt = 0,
			int i_runge_kutta_order = 1,
			double i_simulation_time = -1
	)
	{
		setupBuffers(io_h, i_runge_kutta_order);

		double &dt = o_dt;
		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h,	// input
					io_u,
					io_v,
					*RK_h_t[0],	// output
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_fixed_dt,
					i_simulation_time
			);

			io_h += dt**RK_h_t[0];
			io_u += dt**RK_u_t[0];
			io_v += dt**RK_v_t[0];

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
					io_h,
					io_u,
					io_v,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_fixed_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h + ( dt*a2[0]*(*RK_h_t[0]) ),
					io_u + ( dt*a2[0]*(*RK_u_t[0]) ),
					io_v + ( dt*a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			io_h += dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_h_t[1]) );
			io_u += dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_u_t[1]) );
			io_v += dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_v_t[1]) );
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
					io_h,
					io_u,
					io_v,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_fixed_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( a2[0]*(*RK_h_t[0]) ),
					io_u	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_v	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( a3[0]*(*RK_h_t[0]) + a3[1]*(*RK_h_t[1]) ),
					io_u	+ dt*( a3[0]*(*RK_u_t[0]) + a3[1]*(*RK_u_t[1]) ),
					io_v	+ dt*( a3[0]*(*RK_v_t[0]) + a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			io_h += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) );
			io_u += dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) );
			io_v += dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) );
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
					io_h,
					io_u,
					io_v,
					*RK_h_t[0],
					*RK_u_t[0],
					*RK_v_t[0],
					dt,
					i_fixed_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( a2[0]*(*RK_h_t[0]) ),
					io_u	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_v	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_h_t[1],
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_h_t[1]) ),
					io_u	+ dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_u_t[1]) ),
					io_v	+ dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_v_t[1]) ),
					*RK_h_t[2],
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_t[2]) ),
					io_u	+ dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_t[2]) ),
					io_v	+ dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_t[2]) ),
					*RK_h_t[3],
					*RK_u_t[3],
					*RK_v_t[3],
					dummy_dt,
					dt,
					i_simulation_time + c[2]*dt
			);


			io_h += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) + (b[3]*(*RK_h_t[3])) );
			io_u += dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) + (b[3]*(*RK_u_t[3])) );
			io_v += dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) + (b[3]*(*RK_v_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}
};

#endif
