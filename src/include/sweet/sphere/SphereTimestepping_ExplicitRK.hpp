
#ifndef SPHEREDATA_TIMESTEPPING_EXPLICITRK_HPP__
#define SPHEREDATA_TIMESTEPPING_EXPLICITRK_HPP__

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <limits>

class SphereTimestepping_ExplicitRK
{
	// Runge-Kutta data storages
	std::vector<SphereData_Spectral*> RK_prog0_stage_t;
	std::vector<SphereData_Spectral*> RK_prog1_stage_t;
	std::vector<SphereData_Spectral*> RK_prog2_stage_t;

	std::vector<SphereData_Physical*> RK_vel_u;
	std::vector<SphereData_Physical*> RK_vel_v;

	int runge_kutta_order;

public:
	SphereTimestepping_ExplicitRK()	:
		runge_kutta_order(-1)
	{
	}



	void resetAndSetup(
			const SphereData_Config *i_sphereDataConfig,
			int i_rk_stages			///< Order of Runge-Kutta method
	)
	{
		if (RK_prog0_stage_t.size() != 0)	///< already allocated?
			return;

		runge_kutta_order = i_rk_stages;
		int N = runge_kutta_order;

		if (N <= 0 || N > 4)
			SWEETError("Invalid order for RK time stepping");

		RK_prog0_stage_t.resize(N);
		RK_prog1_stage_t.resize(N);
		RK_prog2_stage_t.resize(N);

		for (int i = 0; i < N; i++)
		{
			RK_prog0_stage_t[i] = new SphereData_Spectral(i_sphereDataConfig);
			RK_prog1_stage_t[i] = new SphereData_Spectral(i_sphereDataConfig);
			RK_prog2_stage_t[i] = new SphereData_Spectral(i_sphereDataConfig);
		}
	}



	void resetAndSetup_na(
			const SphereData_Config *i_sphereDataConfig,
			int i_rk_stages			///< Number of stages
	)
	{
		if (RK_prog0_stage_t.size() != 0)	///< already allocated?
			return;

		runge_kutta_order = i_rk_stages;
		int N = runge_kutta_order;

		if (N <= 0 || N > 4)
			SWEETError("Invalid order for RK time stepping");

		RK_prog0_stage_t.resize(N);
		RK_vel_u.resize(N);
		RK_vel_v.resize(N);

		for (int i = 0; i < N; i++)
		{
			RK_prog0_stage_t[i] = new SphereData_Spectral(i_sphereDataConfig);
			RK_vel_u[i] = new SphereData_Physical(i_sphereDataConfig);
			RK_vel_v[i] = new SphereData_Physical(i_sphereDataConfig);
		}
	}



	~SphereTimestepping_ExplicitRK()
	{
		int N = runge_kutta_order;

		if (RK_prog0_stage_t.size() != 0)
		{
			for (int i = 0; i < N; i++)
				delete RK_prog0_stage_t[i];
			RK_prog0_stage_t.resize(0);
		}

		if (RK_prog1_stage_t.size() != 0)
		{
			for (int i = 0; i < N; i++)
				delete RK_prog1_stage_t[i];
			RK_prog1_stage_t.resize(0);
		}

		if (RK_prog2_stage_t.size() != 0)
		{
			for (int i = 0; i < N; i++)
				delete RK_prog2_stage_t[i];
			RK_prog2_stage_t.resize(0);
		}


		if (RK_vel_u.size() != 0)
		{
			for (int i = 0; i < N; i++)
				delete RK_vel_u[i];
			RK_vel_u.resize(0);
		}

		if (RK_vel_v.size() != 0)
		{
			for (int i = 0; i < N; i++)
				delete RK_vel_v[i];
			RK_vel_v.resize(0);
		}
	}



	/**
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData_Spectral &i_P,	///< prognostic variables
					const SphereData_Spectral &i_u,	///< prognostic variables
					const SphereData_Spectral &i_v,	///< prognostic variables

					SphereData_Spectral &o_P_t,		///< time updates
					SphereData_Spectral &o_u_t,		///< time updates
					SphereData_Spectral &o_v_t,		///< time updates

					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData_Spectral &io_h,
			SphereData_Spectral &io_u,
			SphereData_Spectral &io_v,

			double i_dt,	///< If this value is not equal to 0,
					///< Use this time step size instead of computing one
					///< This also sets o_dt = i_dt

			int i_runge_kutta_order,	///< Order of RK time stepping

			double i_simulation_time	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		resetAndSetup(io_h.sphereDataConfig, i_runge_kutta_order);

		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h,	// input
					io_u,
					io_v,
					*RK_prog0_stage_t[0],	// output
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			io_h += i_dt**RK_prog0_stage_t[0];
			io_u += i_dt**RK_prog1_stage_t[0];
			io_v += i_dt**RK_prog2_stage_t[0];
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h + ( i_dt*a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u + ( i_dt*a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v + ( i_dt*a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			io_h += i_dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_prog0_stage_t[1]) );
			io_u += i_dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_prog1_stage_t[1]) );
			io_v += i_dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_prog2_stage_t[1]) );
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u	+ i_dt*( a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v	+ i_dt*( a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a3[0]*(*RK_prog0_stage_t[0]) + a3[1]*(*RK_prog0_stage_t[1]) ),
					io_u	+ i_dt*( a3[0]*(*RK_prog1_stage_t[0]) + a3[1]*(*RK_prog1_stage_t[1]) ),
					io_v	+ i_dt*( a3[0]*(*RK_prog2_stage_t[0]) + a3[1]*(*RK_prog2_stage_t[1]) ),
					*RK_prog0_stage_t[2],
					*RK_prog1_stage_t[2],
					*RK_prog2_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) );
			io_u += i_dt*( (b[0]*(*RK_prog1_stage_t[0])) + (b[1]*(*RK_prog1_stage_t[1]))  + (b[2]*(*RK_prog1_stage_t[2])) );
			io_v += i_dt*( (b[0]*(*RK_prog2_stage_t[0])) + (b[1]*(*RK_prog2_stage_t[1]))  + (b[2]*(*RK_prog2_stage_t[2])) );
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u	+ i_dt*( a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v	+ i_dt*( a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_prog0_stage_t[1]) ),
					io_u	+ i_dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_prog1_stage_t[1]) ),
					io_v	+ i_dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_prog2_stage_t[1]) ),
					*RK_prog0_stage_t[2],
					*RK_prog1_stage_t[2],
					*RK_prog2_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_prog0_stage_t[2]) ),
					io_u	+ i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_prog1_stage_t[2]) ),
					io_v	+ i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_prog2_stage_t[2]) ),
					*RK_prog0_stage_t[3],
					*RK_prog1_stage_t[3],
					*RK_prog2_stage_t[3],
					i_simulation_time + c[2]*i_dt
			);


			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) + (b[3]*(*RK_prog0_stage_t[3])) );
			io_u += i_dt*( (b[0]*(*RK_prog1_stage_t[0])) + (b[1]*(*RK_prog1_stage_t[1]))  + (b[2]*(*RK_prog1_stage_t[2])) + (b[3]*(*RK_prog1_stage_t[3])) );
			io_v += i_dt*( (b[0]*(*RK_prog2_stage_t[0])) + (b[1]*(*RK_prog2_stage_t[1]))  + (b[2]*(*RK_prog2_stage_t[2])) + (b[3]*(*RK_prog2_stage_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}




	/**
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 *
	 * Special version for Nonlinear Advection term
	 */
	template <class BaseClass>
	void run_timestep_na(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData_Spectral &i_P,	///< prognostic variables
					SphereData_Physical &io_u,
					SphereData_Physical &io_v,

					SphereData_Spectral &o_P_t,		///< time updates

					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData_Spectral &io_h,
			SphereData_Physical &io_u,
			SphereData_Physical &io_v,

			double i_dt,	///< If this value is not equal to 0,
					///< Use this time step size instead of computing one
					///< This also sets o_dt = i_dt

			int i_runge_kutta_order,	///< Order of RK time stepping

			double i_simulation_time	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		resetAndSetup_na(io_h.sphereDataConfig, i_runge_kutta_order);

		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h,	// input
					io_u,
					io_v,
					*RK_prog0_stage_t[0],	// output
					i_simulation_time
			);

			io_h += i_dt**RK_prog0_stage_t[0];
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h + ( i_dt*a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			io_h += i_dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_prog0_stage_t[1]) );
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a3[0]*(*RK_prog0_stage_t[0]) + a3[1]*(*RK_prog0_stage_t[1]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) );
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
					io_h,
					io_u,
					io_v,
					*RK_prog0_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_prog0_stage_t[1]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_prog0_stage_t[2]) ),
					io_u,
					io_v,
					*RK_prog0_stage_t[3],
					i_simulation_time + c[2]*i_dt
			);


			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) + (b[3]*(*RK_prog0_stage_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}




	/**
	 * execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 * This routine is used for the Burgers equation.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData_Spectral &i_u,	///< prognostic variables
					const SphereData_Spectral &i_v,	///< prognostic variables

					SphereData_Spectral &o_u_t,	///< time updates
					SphereData_Spectral &o_v_t,	///< time updates

					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData_Spectral &io_u,
			SphereData_Spectral &io_v,

			double i_dt = 0,	///< If this value is not equal to 0,
						///< Use this time step size instead of computing one
						///< This also sets o_dt = i_dt

			int i_runge_kutta_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u,	// input
					io_v,
					*RK_prog1_stage_t[0],	// output
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			io_u += i_dt**RK_prog1_stage_t[0];
			io_v += i_dt**RK_prog2_stage_t[0];

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
					io_u,
					io_v,
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u + ( i_dt*a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v + ( i_dt*a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			io_u += i_dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_prog1_stage_t[1]) );
			io_v += i_dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_prog2_stage_t[1]) );
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
					io_u,
					io_v,
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ i_dt*( a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v	+ i_dt*( a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ i_dt*( a3[0]*(*RK_prog1_stage_t[0]) + a3[1]*(*RK_prog1_stage_t[1]) ),
					io_v	+ i_dt*( a3[0]*(*RK_prog2_stage_t[0]) + a3[1]*(*RK_prog2_stage_t[1]) ),
					*RK_prog1_stage_t[2],
					*RK_prog2_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			io_u += i_dt*( (b[0]*(*RK_prog1_stage_t[0])) + (b[1]*(*RK_prog1_stage_t[1]))  + (b[2]*(*RK_prog1_stage_t[2])) );
			io_v += i_dt*( (b[0]*(*RK_prog2_stage_t[0])) + (b[1]*(*RK_prog2_stage_t[1]))  + (b[2]*(*RK_prog2_stage_t[2])) );
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
					io_u,
					io_v,
					*RK_prog1_stage_t[0],
					*RK_prog2_stage_t[0],
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ i_dt*( a2[0]*(*RK_prog1_stage_t[0]) ),
					io_v	+ i_dt*( a2[0]*(*RK_prog2_stage_t[0]) ),
					*RK_prog1_stage_t[1],
					*RK_prog2_stage_t[1],
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ i_dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_prog1_stage_t[1]) ),
					io_v	+ i_dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_prog2_stage_t[1]) ),
					*RK_prog1_stage_t[2],
					*RK_prog2_stage_t[2],
					i_simulation_time + c[1]*i_dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ i_dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_prog1_stage_t[2]) ),
					io_v	+ i_dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_prog2_stage_t[2]) ),
					*RK_prog1_stage_t[3],
					*RK_prog2_stage_t[3],
					i_simulation_time + c[2]*i_dt
			);


			io_u += i_dt*( (b[0]*(*RK_prog1_stage_t[0])) + (b[1]*(*RK_prog1_stage_t[1]))  + (b[2]*(*RK_prog1_stage_t[2])) + (b[3]*(*RK_prog1_stage_t[3])) );
			io_v += i_dt*( (b[0]*(*RK_prog2_stage_t[0])) + (b[1]*(*RK_prog2_stage_t[1]))  + (b[2]*(*RK_prog2_stage_t[2])) + (b[3]*(*RK_prog2_stage_t[3])) );
		}
		else
		{
			SWEETError("This order of the Runge-Kutta time stepping is not supported!");
		}
	}



	/**
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData_Spectral &i_h,		///< prognostic variables
					SphereData_Spectral &o_h_t,			///< time updates

					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData_Spectral &io_h,

			double i_dt,				///< If this value is not equal to 0,
										///< Use this time step size instead of computing one
										///< This also sets o_dt = i_dt

			int i_runge_kutta_order,	///< Order of RK time stepping

			double i_simulation_time	///< Current simulation time.
										///< This gets e.g. important for tidal waves
	)
	{
//		resetAndSetup(io_h, i_runge_kutta_order);

		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h,	// input
					*RK_prog0_stage_t[0],	// output
					i_dt,
					i_simulation_time
			);

			io_h += i_dt**RK_prog0_stage_t[0];

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
					io_h,
					*RK_prog0_stage_t[0],
					i_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h + ( i_dt*a2[0]*(*RK_prog0_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					i_dt,
					i_simulation_time + c[0]*i_dt
			);

			io_h += i_dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_prog0_stage_t[1]) );
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
					io_h,
					*RK_prog0_stage_t[0],
					i_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					i_dt,
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a3[0]*(*RK_prog0_stage_t[0]) + a3[1]*(*RK_prog0_stage_t[1]) ),
					*RK_prog0_stage_t[2],
					i_dt,
					i_simulation_time + c[1]*i_dt
			);

			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) );
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
					*RK_prog0_stage_t[0],
					i_dt,
					i_simulation_time
			);

			// STAGE 2
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( a2[0]*(*RK_prog0_stage_t[0]) ),
					*RK_prog0_stage_t[1],
					i_dt,
					i_simulation_time + c[0]*i_dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_prog0_stage_t[1]) ),
					*RK_prog0_stage_t[2],
					i_dt,
					i_simulation_time + c[1]*i_dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ i_dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_prog0_stage_t[2]) ),
					*RK_prog0_stage_t[3],
					i_dt,
					i_simulation_time + c[2]*i_dt
			);


			io_h += i_dt*( (b[0]*(*RK_prog0_stage_t[0])) + (b[1]*(*RK_prog0_stage_t[1]))  + (b[2]*(*RK_prog0_stage_t[2])) + (b[3]*(*RK_prog0_stage_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}


};

#endif
