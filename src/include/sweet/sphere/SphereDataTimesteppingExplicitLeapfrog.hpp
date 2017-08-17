
#ifndef TIMESTEPPING_EXPLICIT_LEAPFROG_HPP
#define TIMESTEPPING_EXPLICIT_LEAPFROG_HPP

#include <sweet/sphere/SphereData.hpp>
#include <limits>



class SphereDataTimesteppingExplicitLeapfrog
{
	SphereDataConfig *sphereDataConfig;

	// Previous time step values
	SphereData RK_h_prev;
	SphereData RK_u_prev;
	SphereData RK_v_prev;

	// time step tendencies
	SphereData RK_h_dt;
	SphereData RK_u_dt;
	SphereData RK_v_dt;

	// Temporary time step values
	SphereData RK_h_tmp;
	SphereData RK_u_tmp;
	SphereData RK_v_tmp;


	/**
	 * Robert-Asselin filter coefficient
	 */
	double leapfrog_robert_asselin_filter;

	int leapfrog_order;
	int timestep_id;

public:
	SphereDataTimesteppingExplicitLeapfrog(
			SphereDataConfig *i_sphereDataconfig
	)	:
		sphereDataConfig(i_sphereDataconfig),

		RK_h_prev(i_sphereDataconfig),
		RK_u_prev(i_sphereDataconfig),
		RK_v_prev(i_sphereDataconfig),

		RK_h_dt(i_sphereDataconfig),
		RK_u_dt(i_sphereDataconfig),
		RK_v_dt(i_sphereDataconfig),

		RK_h_tmp(i_sphereDataconfig),
		RK_u_tmp(i_sphereDataconfig),
		RK_v_tmp(i_sphereDataconfig),

		leapfrog_robert_asselin_filter(0),
		timestep_id(0)
	{
	}



	void setup(
			double i_leapfrog_robert_asselin_filter = 0	/// Filter value for Robert Asselin filter. Value of 0 means no filter
	)
	{
		// reset the time step id
		timestep_id = 0;
		leapfrog_robert_asselin_filter = i_leapfrog_robert_asselin_filter;
	}



	~SphereDataTimesteppingExplicitLeapfrog()
	{
	}



	/**
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData &i_P,	///< prognostic variables
					const SphereData &i_u,	///< prognostic variables
					const SphereData &i_v,	///< prognostic variables

					SphereData &o_P_t,		///< time updates
					SphereData &o_u_t,		///< time updates
					SphereData &o_v_t,		///< time updates

					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_h,
			SphereData &io_u,
			SphereData &io_v,

			double i_use_fixed_dt = 0,		///< If this value is not equal to 0,
											///< Use this time step size instead of computing one
											///< This also sets o_dt = i_use_fixed_dt

			int i_runge_kutta_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		double &dt = i_use_fixed_dt;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_h,		// input
				io_u,
				io_v,
				RK_h_dt,	// output
				RK_u_dt,
				RK_v_dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				RK_h_prev = io_h;
				RK_u_prev = io_u;
				RK_v_prev = io_v;

#if 0
				// do standard Euler time step
				io_h += dt * RK_h_dt;
				io_u += dt * RK_u_dt;
				io_v += dt * RK_v_dt;
#else

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

				SphereData h_t0(io_h.sphereDataConfig);
				SphereData u_t0(io_h.sphereDataConfig);
				SphereData v_t0(io_h.sphereDataConfig);

				// STAGE 1
				(i_baseClass->*i_compute_euler_timestep_update)(
						io_h,
						io_u,
						io_v,
						h_t0,//*RK_h_t[0],
						u_t0,//*RK_u_t[0],
						v_t0,//*RK_v_t[0],
						i_use_fixed_dt,
						i_simulation_time
				);

				SphereData h_t1(io_h.sphereDataConfig);
				SphereData u_t1(io_h.sphereDataConfig);
				SphereData v_t1(io_h.sphereDataConfig);

				// STAGE 2
				(i_baseClass->*i_compute_euler_timestep_update)(
						io_h + ( dt*a2[0]*h_t0),//(*RK_h_t[0]) ),
						io_u + ( dt*a2[0]*u_t0),//(*RK_u_t[0]) ),
						io_v + ( dt*a2[0]*v_t0),//(*RK_v_t[0]) ),
						h_t1,//*RK_h_t[1],
						u_t1,//*RK_u_t[1],
						v_t1,//*RK_v_t[1],
						dt,
						i_simulation_time + c[0]*dt
				);

				//io_h += dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*(*RK_h_t[1]) );
				//io_u += dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*(*RK_u_t[1]) );
				//io_v += dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*(*RK_v_t[1]) );
				io_h += dt*(/* b[0]*(*RK_h_t[0]) +*/ b[1]*h_t1 );
				io_u += dt*(/* b[0]*(*RK_u_t[0]) +*/ b[1]*u_t1 );
				io_v += dt*(/* b[0]*(*RK_v_t[0]) +*/ b[1]*v_t1 );
#endif
			}
			else
			{
				RK_h_tmp.swap(io_h);
				RK_u_tmp.swap(io_u);
				RK_v_tmp.swap(io_v);

				io_h = RK_h_prev + 2.0*dt * RK_h_dt;
				io_u = RK_u_prev + 2.0*dt * RK_u_dt;
				io_v = RK_v_prev + 2.0*dt * RK_v_dt;

				RK_h_prev.swap(RK_h_tmp);
				RK_u_prev.swap(RK_u_tmp);
				RK_v_prev.swap(RK_v_tmp);
			}
		}
		else
		{
			/*
			 * See Robert(1966)
			 *
			 * F*(t+dt) = F(t-dt) + 2 dt (dF/dt)*
			 *
			 * F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
			 *
			 * "*" denote preliminary values
			 */
			if (timestep_id == 0)
			{
				// backup previous values
				RK_h_prev = io_h;
				RK_u_prev = io_u;
				RK_v_prev = io_v;

				// do standard Euler time step
				io_h += dt * RK_h_dt;
				io_u += dt * RK_u_dt;
				io_v += dt * RK_v_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)

				// V(t+dt) = U(t-dt) + 2 dt (dU/dt)
				RK_h_tmp = RK_h_prev + 2.0*dt*RK_h_dt;
				RK_u_tmp = RK_u_prev + 2.0*dt*RK_u_dt;
				RK_v_tmp = RK_v_prev + 2.0*dt*RK_v_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]

				// U(t) = V(t) + 0.5*v*[ V(t+dt) + U(t-dt) - 2 V(t) ]
				RK_h_prev = io_h + leapfrog_robert_asselin_filter*0.5*(RK_h_tmp + RK_h_prev - 2.0*io_h);
				RK_u_prev = io_u + leapfrog_robert_asselin_filter*0.5*(RK_u_tmp + RK_u_prev - 2.0*io_u);
				RK_v_prev = io_v + leapfrog_robert_asselin_filter*0.5*(RK_v_tmp + RK_v_prev - 2.0*io_v);

				io_h.swap(RK_h_tmp);
				io_u.swap(RK_u_tmp);
				io_v.swap(RK_v_tmp);
			}
		}

		timestep_id++;
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
					const SphereData &i_u,	///< prognostic variables
					const SphereData &i_v,	///< prognostic variables

					SphereData &o_u_t,	///< time updates
					SphereData &o_v_t,	///< time updates

					double i_use_fixed_dt,	///< if this value is not equal to 0,
											///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_u,
			SphereData &io_v,

			double i_use_fixed_dt = 0,		///< If this value is not equal to 0,
											///< Use this time step size instead of computing one
											///< This also sets o_dt = i_use_fixed_dt

			int i_leapfrog_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		/*
		 * See
		 * "Analysis of time filters used with the leapfrog scheme"
		 * Yong Li, Catalin Trenchea
		 */

		double &dt = i_use_fixed_dt;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_u,		// input
				io_v,		// input
				RK_u_dt,	// output
				RK_v_dt,	// output
				dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				RK_u_prev = io_u;
				RK_v_prev = io_v;

				// do standard Euler time step
				io_u += dt * RK_u_dt;
				io_v += dt * RK_v_dt;
			}
			else
			{
				RK_u_tmp.swap(io_u);
				RK_v_tmp.swap(io_v);

				io_u = RK_u_prev + 2.0*dt * RK_u_dt;
				io_v = RK_v_prev + 2.0*dt * RK_v_dt;

				RK_u_prev.swap(RK_u_tmp);
				RK_v_prev.swap(RK_v_tmp);
			}
		}
		else
		{
			/*
			 * See Robert(1966)
			 *
			 * F*(t+dt) = F(t-dt) + 2 dt (dF/dt)*
			 *
			 * F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
			 *
			 * "*" denote preliminary values
			 */
			if (timestep_id == 0)
			{
				// backup previous values
				RK_u_prev = io_u;
				RK_v_prev = io_v;

				// do standard Euler time step
				io_u += dt * RK_u_dt;
				io_v += dt * RK_v_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)
				RK_u_tmp = RK_u_prev + 2.0*dt*RK_u_dt;
				RK_v_tmp = RK_v_prev + 2.0*dt*RK_v_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
				RK_u_prev = io_u + leapfrog_robert_asselin_filter*0.5*(RK_u_tmp + RK_u_prev - 2.0*io_u);
				RK_v_prev = io_v + leapfrog_robert_asselin_filter*0.5*(RK_v_tmp + RK_v_prev - 2.0*io_v);

				io_u.swap(RK_u_tmp);
				io_v.swap(RK_v_tmp);
			}
		}

		timestep_id++;
	}



	/**
	 * Execute a Runge-Kutta timestep with the order
	 * specified in the simulation variables.
	 */
	template <class BaseClass>
	void run_timestep(
			BaseClass *i_baseClass,
			void (BaseClass::*i_compute_euler_timestep_update)(
					const SphereData &i_h,		///< prognostic variables
					SphereData &o_h_t,			///< time updates

					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_h,

			double i_use_fixed_dt = 0,		///< If this value is not equal to 0,
											///< Use this time step size instead of computing one
											///< This also sets o_dt = i_use_fixed_dt

			int i_runge_kutta_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1	///< Current simulation time.
											///< This gets e.g. important for tidal waves
	)
	{
		double &dt = i_use_fixed_dt;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_h,		// input
				RK_h_dt,	// output
				dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				RK_h_prev = io_h;

				// do standard Euler time step
				io_h += dt * RK_h_dt;
			}
			else
			{
				RK_h_tmp.swap(io_h);

				io_h = RK_h_prev + 2.0*dt * RK_h_dt;

				RK_h_prev.swap(RK_h_tmp);
			}
		}
		else
		{
			if (timestep_id == 0)
			{
				// backup previous values
				RK_h_prev = io_h;

				// do standard Euler time step
				io_h += dt * RK_h_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)
				RK_h_tmp = RK_h_prev + 2.0*dt*RK_h_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
				RK_h_prev = io_h + leapfrog_robert_asselin_filter*0.5*(RK_h_tmp + RK_h_prev - 2.0*io_h);

				io_h.swap(RK_h_tmp);
			}
		}

		timestep_id++;
	}

};

#endif
