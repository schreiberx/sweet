
#ifndef SPHEREDATAPHYSICAL_TIMESTEPPING_EXPLICIT_RK_HPP
#define SPHEREDATAPHYSICAL_TIMESTEPPING_EXPLICIT_RK_HPP

#include <sweet/sphere/SphereDataPhysical.hpp>
#include <limits>

class SphereDataPhysicalTimesteppingExplicitRK
{
	// Runge-Kutta data storages
	SphereDataPhysical** RK_h_t;
	SphereDataPhysical** RK_u_t;
	SphereDataPhysical** RK_v_t;

	int runge_kutta_order;


public:
	SphereDataPhysicalTimesteppingExplicitRK()	:
		RK_h_t(nullptr),
		RK_u_t(nullptr),
		RK_v_t(nullptr),
		runge_kutta_order(-1)
	{
	}



	void resetAndSetup(
			SphereDataConfig *i_sphereDataConfig,
			int i_rk_order			///< Order of Runge-Kutta method
	)
	{
		if (RK_h_t != nullptr)	///< already allocated?
			return;

		runge_kutta_order = i_rk_order;
		int N = i_rk_order;

		if (N <= 0 || N > 4)
			FatalError("Invalid order for RK time stepping");

		RK_h_t = new SphereDataPhysical*[N];
		RK_u_t = new SphereDataPhysical*[N];
		RK_v_t = new SphereDataPhysical*[N];

		for (int i = 0; i < N; i++)
		{
			RK_h_t[i] = new SphereDataPhysical(i_sphereDataConfig);
			RK_u_t[i] = new SphereDataPhysical(i_sphereDataConfig);
			RK_v_t[i] = new SphereDataPhysical(i_sphereDataConfig);
		}
	}



	~SphereDataPhysicalTimesteppingExplicitRK()
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
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereDataPhysical &i_P,	///< prognostic variables
					const SphereDataPhysical &i_u,	///< prognostic variables
					const SphereDataPhysical &i_v,	///< prognostic variables

					SphereDataPhysical &o_P_t,		///< time updates
					SphereDataPhysical &o_u_t,		///< time updates
					SphereDataPhysical &o_v_t,		///< time updates

					double &o_dt,				///< time step restriction
					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereDataPhysical &io_h,
			SphereDataPhysical &io_u,
			SphereDataPhysical &io_v,

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
//		resetAndSetup(io_h, i_runge_kutta_order);

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
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;


#if 0
			io_h = io_h.getSphereDataPhysicalPhysical() + dt*(*RK_h_t[0]).getSphereDataPhysicalPhysical();
			io_u = io_u.getSphereDataPhysicalPhysical() + dt*(*RK_u_t[0]).getSphereDataPhysicalPhysical();
			io_v = io_v.getSphereDataPhysicalPhysical() + dt*(*RK_v_t[0]).getSphereDataPhysicalPhysical();
#else
			io_h += dt**RK_h_t[0];
			io_u += dt**RK_u_t[0];
			io_v += dt**RK_v_t[0];
#endif
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
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

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
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

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
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

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




	/**
	 * execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 * This routine is used for the Burgers equation.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereDataPhysical &i_u,	///< prognostic variables
					const SphereDataPhysical &i_v,	///< prognostic variables

					SphereDataPhysical &o_u_t,	///< time updates
					SphereDataPhysical &o_v_t,	///< time updates

					double &o_dt,			///< time step restriction
					double i_use_fixed_dt,	///< if this value is not equal to 0,
											///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereDataPhysical &io_u,
			SphereDataPhysical &io_v,

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
//		resetAndSetup(io_u, i_runge_kutta_order);

		double &dt = o_dt;
		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u,	// input
					io_v,
					*RK_u_t[0],	// output
					*RK_v_t[0],
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;

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
					io_u,
					io_v,
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
					io_u + ( dt*a2[0]*(*RK_u_t[0]) ),
					io_v + ( dt*a2[0]*(*RK_v_t[0]) ),
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

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
					io_u,
					io_v,
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
					io_u	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_v	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ dt*( a3[0]*(*RK_u_t[0]) + a3[1]*(*RK_u_t[1]) ),
					io_v	+ dt*( a3[0]*(*RK_v_t[0]) + a3[1]*(*RK_v_t[1]) ),
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

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
					io_u,
					io_v,
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
					io_u	+ dt*( a2[0]*(*RK_u_t[0]) ),
					io_v	+ dt*( a2[0]*(*RK_v_t[0]) ),
					*RK_u_t[1],
					*RK_v_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ dt*( /*a3[0]*(*RK_u_t[0]) +*/ a3[1]*(*RK_u_t[1]) ),
					io_v	+ dt*( /*a3[0]*(*RK_v_t[0]) +*/ a3[1]*(*RK_v_t[1]) ),
					*RK_u_t[2],
					*RK_v_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_u	+ dt*( /*a4[0]*(*RK_u_t[0]) + a4[1]*(*RK_u_t[1]) +*/ a4[2]*(*RK_u_t[2]) ),
					io_v	+ dt*( /*a4[0]*(*RK_v_t[0]) + a4[1]*(*RK_v_t[1]) +*/ a4[2]*(*RK_v_t[2]) ),
					*RK_u_t[3],
					*RK_v_t[3],
					dummy_dt,
					dt,
					i_simulation_time + c[2]*dt
			);


			io_u += dt*( (b[0]*(*RK_u_t[0])) + (b[1]*(*RK_u_t[1]))  + (b[2]*(*RK_u_t[2])) + (b[3]*(*RK_u_t[3])) );
			io_v += dt*( (b[0]*(*RK_v_t[0])) + (b[1]*(*RK_v_t[1]))  + (b[2]*(*RK_v_t[2])) + (b[3]*(*RK_v_t[3])) );
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
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereDataPhysical &i_h,		///< prognostic variables
					SphereDataPhysical &o_h_t,			///< time updates

					double &o_dt,				///< time step restriction
					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereDataPhysical &io_h,

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
//		resetAndSetup(io_h, i_runge_kutta_order);

		double &dt = o_dt;
		if (i_runge_kutta_order == 1)
		{
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h,	// input
					*RK_h_t[0],	// output
					dt,
					i_use_fixed_dt,
					i_simulation_time
			);

			// padding to max simulation time if exceeding the maximum
			if (i_max_simulation_time >= 0)
				if (dt+i_simulation_time > i_max_simulation_time)
					dt = i_max_simulation_time-i_simulation_time;


			io_h += dt**RK_h_t[0];

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
					*RK_h_t[0],
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
					io_h + ( dt*a2[0]*(*RK_h_t[0]) ),
					*RK_h_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			io_h += dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_h_t[1]) );
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
					*RK_h_t[0],
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
					io_h	+ dt*( a2[0]*(*RK_h_t[0]) ),
					*RK_h_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( a3[0]*(*RK_h_t[0]) + a3[1]*(*RK_h_t[1]) ),
					*RK_h_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			io_h += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) );
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
					*RK_h_t[0],
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
					io_h	+ dt*( a2[0]*(*RK_h_t[0]) ),
					*RK_h_t[1],
					dummy_dt,
					dt,
					i_simulation_time + c[0]*dt
			);

			// STAGE 3
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( /*a3[0]*(*RK_P_t[0]) +*/ a3[1]*(*RK_h_t[1]) ),
					*RK_h_t[2],
					dummy_dt,
					dt,
					i_simulation_time + c[1]*dt
			);

			// STAGE 4
			(i_baseClass->*i_compute_euler_timestep_update)(
					io_h	+ dt*( /*a4[0]*(*RK_P_t[0]) + a4[1]*(*RK_P_t[1]) +*/ a4[2]*(*RK_h_t[2]) ),
					*RK_h_t[3],
					dummy_dt,
					dt,
					i_simulation_time + c[2]*dt
			);


			io_h += dt*( (b[0]*(*RK_h_t[0])) + (b[1]*(*RK_h_t[1]))  + (b[2]*(*RK_h_t[2])) + (b[3]*(*RK_h_t[3])) );
		}
		else
		{
			std::cerr << "This order of the Runge-Kutta time stepping is not supported!" << std::endl;
			exit(-1);
		}
	}


};

#endif
