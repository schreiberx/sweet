
#ifndef TIMESTEPPING_EXPLICIT_LEAPFROG_HPP
#define TIMESTEPPING_EXPLICIT_LEAPFROG_HPP

#include <sweet/sphere/SphereData.hpp>
#include <limits>



class SphereDataTimesteppingExplicitLeapfrog
{
	// Previous time step values
	SphereData* RK_h_prev;
	SphereData* RK_u_prev;
	SphereData* RK_v_prev;

	// time step tendencies
	SphereData* RK_h_dt;
	SphereData* RK_u_dt;
	SphereData* RK_v_dt;

	// Temporary time step values
	SphereData* RK_h_tmp;
	SphereData* RK_u_tmp;
	SphereData* RK_v_tmp;


	/**
	 * Robert-Asselin filter coefficient
	 */
	double leapfrog_robert_asselin_filter;

	int leapfrog_order;
	int timestep_id;

public:
	SphereDataTimesteppingExplicitLeapfrog()	:
		RK_h_prev(nullptr),
		RK_u_prev(nullptr),
		RK_v_prev(nullptr),
		leapfrog_robert_asselin_filter(0),
		leapfrog_order(-1),
		timestep_id(0)
	{
	}


	void cleanup()
	{
		if (RK_h_prev)
		{
			delete RK_h_prev;
			delete RK_u_prev;
			delete RK_v_prev;

			RK_h_prev = nullptr;
			RK_u_prev = nullptr;
			RK_v_prev = nullptr;

			delete RK_h_dt;
			delete RK_u_dt;
			delete RK_v_dt;

			RK_h_dt = nullptr;
			RK_u_dt = nullptr;
			RK_v_dt = nullptr;

			delete RK_h_tmp;
			delete RK_u_tmp;
			delete RK_v_tmp;

			RK_h_tmp = nullptr;
			RK_u_tmp = nullptr;
			RK_v_tmp = nullptr;
		}
	}


	void resetAndSetup(
			const SphereData &i_test_buffer,	///< array of example data to know dimensions of buffers
			int i_leapfrog_order,			///< Order of Leapfrog method
			double i_leapfrog_robert_asselin_filter = 0	/// Filter value for Robert Asselin filter. Value of 0 means no filter
	)
	{
		// reset the time step id
		timestep_id = 0;
		leapfrog_order = i_leapfrog_order;
		leapfrog_robert_asselin_filter = i_leapfrog_robert_asselin_filter;

		if (RK_h_prev != nullptr)	///< already allocated?
			return;

		int N = i_leapfrog_order;

		if (N <= 0 || N > 1)
			FatalError("Only 1st order leapfrog is currently supported!");


		// storage for previous time step
		RK_h_prev = new SphereData(i_test_buffer.sphereDataConfig);
		RK_u_prev = new SphereData(i_test_buffer.sphereDataConfig);
		RK_v_prev = new SphereData(i_test_buffer.sphereDataConfig);

		// storage for time step tendencies
		RK_h_dt = new SphereData(i_test_buffer.sphereDataConfig);
		RK_u_dt = new SphereData(i_test_buffer.sphereDataConfig);
		RK_v_dt = new SphereData(i_test_buffer.sphereDataConfig);

		// storage for temporary variables
		RK_h_tmp = new SphereData(i_test_buffer.sphereDataConfig);
		RK_u_tmp = new SphereData(i_test_buffer.sphereDataConfig);
		RK_v_tmp = new SphereData(i_test_buffer.sphereDataConfig);
	}



	~SphereDataTimesteppingExplicitLeapfrog()
	{
		if (RK_h_prev != nullptr)
		{
			cleanup();
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
					const SphereData &i_P,	///< prognostic variables
					const SphereData &i_u,	///< prognostic variables
					const SphereData &i_v,	///< prognostic variables

					SphereData &o_P_t,		///< time updates
					SphereData &o_u_t,		///< time updates
					SphereData &o_v_t,		///< time updates

					double &o_dt,				///< time step restriction
					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_h,
			SphereData &io_u,
			SphereData &io_v,

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


		double &dt = o_dt;

		SphereData &h_prev_u = *RK_h_prev;
		SphereData &u_prev_u = *RK_u_prev;
		SphereData &v_prev_u = *RK_v_prev;

		SphereData &h_dt = *RK_h_dt;
		SphereData &u_dt = *RK_u_dt;
		SphereData &v_dt = *RK_v_dt;

		SphereData &h_tmp = *RK_h_tmp;
		SphereData &u_tmp = *RK_u_tmp;
		SphereData &v_tmp = *RK_v_tmp;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_h,		// input
				io_u,
				io_v,
				h_dt,	// output
				u_dt,
				v_dt,
				dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				h_prev_u = io_h;
				u_prev_u = io_u;
				v_prev_u = io_v;

				// do standard Euler time step
				io_h += dt * h_dt;
				io_u += dt * u_dt;
				io_v += dt * v_dt;

//				std::cout << SphereData(io_h).physical_reduce_sum() << "\t" << SphereData(io_u).physical_reduce_sum() << "\t" << SphereData(io_v).physical_reduce_sum() << std::endl;
			}
			else
			{
				h_tmp.swap(io_h);
				u_tmp.swap(io_u);
				v_tmp.swap(io_v);

				io_h = h_prev_u + 2.0*dt * h_dt;
				io_u = u_prev_u + 2.0*dt * u_dt;
				io_v = v_prev_u + 2.0*dt * v_dt;

				h_prev_u.swap(h_tmp);
				u_prev_u.swap(u_tmp);
				v_prev_u.swap(v_tmp);
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
				h_prev_u = io_h;
				u_prev_u = io_u;
				v_prev_u = io_v;

				// do standard Euler time step
				io_h += dt * h_dt;
				io_u += dt * u_dt;
				io_v += dt * v_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)

				// V(t+dt) = U(t-dt) + 2 dt (dU/dt)
				h_tmp = h_prev_u + 2.0*dt*h_dt;
				u_tmp = u_prev_u + 2.0*dt*u_dt;
				v_tmp = v_prev_u + 2.0*dt*v_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]

				// U(t) = V(t) + 0.5*v*[ V(t+dt) + U(t-dt) - 2 V(t) ]
				h_prev_u = io_h + leapfrog_robert_asselin_filter*0.5*(h_tmp + h_prev_u - 2.0*io_h);
				u_prev_u = io_u + leapfrog_robert_asselin_filter*0.5*(u_tmp + u_prev_u - 2.0*io_u);
				v_prev_u = io_v + leapfrog_robert_asselin_filter*0.5*(v_tmp + v_prev_u - 2.0*io_v);

				io_h.swap(h_tmp);
				io_u.swap(u_tmp);
				io_v.swap(v_tmp);
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

					double &o_dt,			///< time step restriction
					double i_use_fixed_dt,	///< if this value is not equal to 0,
											///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_u,
			SphereData &io_v,

			double &o_dt,					///< return time step size for the computed time step

			double i_use_fixed_dt = 0,		///< If this value is not equal to 0,
											///< Use this time step size instead of computing one
											///< This also sets o_dt = i_use_fixed_dt

			int i_leapfrog_order = 1,	///< Order of RK time stepping

			double i_simulation_time = -1,	///< Current simulation time.
											///< This gets e.g. important for tidal waves

			double i_max_simulation_time = std::numeric_limits<double>::infinity()	///< limit the maximum simulation time
	)
	{
		/*
		 * See
		 * "Analysis of time filters used with the leapfrog scheme"
		 * Yong Li, Catalin Trenchea
		 */


		double &dt = o_dt;

		SphereData &u_prev = *RK_u_prev;
		SphereData &v_prev = *RK_v_prev;

		SphereData &u_dt = *RK_u_dt;
		SphereData &v_dt = *RK_v_dt;

		SphereData &u_tmp = *RK_u_tmp;
		SphereData &v_tmp = *RK_v_tmp;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_u,		// input
				io_v,		// input
				u_dt,	// output
				v_dt,	// output
				dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				u_prev = io_u;
				v_prev = io_v;

				// do standard Euler time step
				io_u += dt * u_dt;
				io_v += dt * v_dt;
			}
			else
			{
				u_tmp.swap(io_u);
				v_tmp.swap(io_v);

				io_u = u_prev + 2.0*dt * u_dt;
				io_v = v_prev + 2.0*dt * v_dt;

				u_prev.swap(u_tmp);
				v_prev.swap(v_tmp);
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
				u_prev = io_u;
				v_prev = io_v;

				// do standard Euler time step
				io_u += dt * u_dt;
				io_v += dt * v_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)
				u_tmp = u_prev + 2.0*dt*u_dt;
				v_tmp = v_prev + 2.0*dt*v_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
				u_prev = io_u + leapfrog_robert_asselin_filter*0.5*(u_tmp + u_prev - 2.0*io_u);
				v_prev = io_v + leapfrog_robert_asselin_filter*0.5*(v_tmp + v_prev - 2.0*io_v);

				io_u.swap(u_tmp);
				io_v.swap(v_tmp);
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

					double &o_dt,				///< time step restriction
					double i_use_fixed_dt,		///< if this value is not equal to 0,
												///< use this time step size instead of computing one
					double i_simulation_time	///< simulation time, e.g. for tidal waves
			),

			SphereData &io_h,

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
		double &dt = o_dt;

		SphereData &h_prev = *RK_h_prev;
		SphereData &h_dt = *RK_h_dt;
		SphereData &h_tmp = *RK_h_tmp;

		(i_baseClass->*i_compute_euler_timestep_update)(
				io_h,		// input
				h_dt,	// output
				dt,
				i_use_fixed_dt,
				i_simulation_time
		);

		if (leapfrog_robert_asselin_filter == 0)
		{
			if (timestep_id == 0)
			{
				h_prev = io_h;

				// do standard Euler time step
				io_h += dt * h_dt;
			}
			else
			{
				h_tmp.swap(io_h);

				io_h = h_prev + 2.0*dt * h_dt;

				h_prev.swap(h_tmp);
			}
		}
		else
		{
			if (timestep_id == 0)
			{
				// backup previous values
				h_prev = io_h;

				// do standard Euler time step
				io_h += dt * h_dt;
			}
			else
			{
				// F*(t+dt) = F(t-dt) + 2 dt (dF/dt)
				h_tmp = h_prev + 2.0*dt*h_dt;

				// F(t) = F*(t) + v*[ f*(t+dt) + F(t-dt) - 2 F*(t) ]
				h_prev = io_h + leapfrog_robert_asselin_filter*0.5*(h_tmp + h_prev - 2.0*io_h);

				io_h.swap(h_tmp);
			}
		}

		timestep_id++;
	}

};

#endif
