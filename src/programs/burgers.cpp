/*
 * Burgers equation
 *
 */

#include "../include/sweet/plane/PlaneData.hpp"
#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <benchmarks_plane/BurgersValidationBenchmarks.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>
#include <sweet/plane/Convert_PlaneData_to_ScalarDataArray.hpp>
#include <sweet/plane/Convert_ScalarDataArray_to_PlaneData.hpp>
#include <sweet/FatalError.hpp>

#include "burgers/burgers_HelmholtzSolver.hpp"

#include <sweet/Stopwatch.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 0
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal.hpp>
#endif


#ifndef SWEET_PARAREAL
#error "ACTIVATE PARAREAL compile option!"
#endif



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

// Input parameters (cmd line)

//general parameters
SimulationVariables simVars;

const int NUM_OF_UNKNOWNS=2;

//specific parameters
bool param_compute_error = 0;
bool param_use_staggering = 0;
bool param_semilagrangian = 0;


class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	PlaneData prog_u, prog_v;

	// Prognostic variables at time step t-dt
	PlaneData prog_u_prev, prog_v_prev;

	// temporary variables - may be overwritten, use locally
	PlaneData tmp;

	// Points mapping [0,simVars.sim.domain_size[0])x[0,simVars.sim.domain_size[1])
	// with resolution simVars.sim.resolution
	ScalarDataArray pos_x, pos_y;

	// Arrival points for semi-lag
	ScalarDataArray posx_a, posy_a;

	// Departure points for semi-lag
	ScalarDataArray posx_d, posy_d;

	//Staggering displacement array (use 0.5 for each displacement)
	// [0] - delta x of u variable
	// [1] - delta y of u variable
	// [2] - delta x of v variable
	// [3] - delta y of v variable
	// Default - A grid (there is shift in x,y of 1/2 to all vars)
	// For C grid use {0,-0.5,-0.5,0}
	double stag_displacement[4] = {-0.5,-0.5,-0.5,-0.5};
	double stag_h[2] = {-0.5,-0.5};
	double stag_u[2] = {-0.5,-0.5};
	double stag_v[2] = {-0.5,-0.5};

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	// Max difference to initial conditions
	double benchmark_diff_u;
	double benchmark_diff_v;

	// Error to analytical solution, if it exists
	PlaneData benchmark_analytical_error;

	// Error measures L2 norm
	double benchmark_analytical_error_rms_u;
	double benchmark_analytical_error_rms_v;

	// Error measures max norm
	double benchmark_analytical_error_maxabs_u;
	double benchmark_analytical_error_maxabs_v;

	// Finite difference operators
	PlaneOperators op;

	// Runge-Kutta stuff
	PlaneDataTimesteppingRK timestepping;

	// Interpolation stuff
	PlaneDataSampler sampler2D;

	// Semi-Lag stuff
	SemiLagrangian semiLagrangian;


	/**
	 * Two dimensional Burgers equation
	 *
	 * Equations:
	 *
	 *     \f$ u_t + uu_x + vu_y = \nu (u_xx + u_yy) \f$
	 *     \f$ v_t + uv_x + vv_y = \nu (v_xx + v_yy) \f$
	 *
	 *   ______________
	 *   |            |
	 *   |    u0,1    |
	 *   v0,0 P0,0 v1,0
	 *   |    u0,0    |
	 *   |____________|
	 */
public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup
		prog_u(planeDataConfig),	// velocity (x-direction)
		prog_v(planeDataConfig),	// velocity (y-direction)

		prog_u_prev(planeDataConfig),
		prog_v_prev(planeDataConfig),

		tmp(planeDataConfig),

		pos_x(planeDataConfig->physical_array_data_number_of_elements),
		pos_y(planeDataConfig->physical_array_data_number_of_elements),

		posx_a(planeDataConfig->physical_array_data_number_of_elements),
		posy_a(planeDataConfig->physical_array_data_number_of_elements),

		posx_d(planeDataConfig->physical_array_data_number_of_elements),
		posy_d(planeDataConfig->physical_array_data_number_of_elements),

		benchmark_analytical_error(planeDataConfig),

		// Init operators
		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_u(planeDataConfig), _parareal_data_start_v(planeDataConfig),
		_parareal_data_start_u_prev(planeDataConfig), _parareal_data_start_v_prev(planeDataConfig),
		_parareal_data_fine_u(planeDataConfig), _parareal_data_fine_v(planeDataConfig),
		_parareal_data_coarse_u(planeDataConfig), _parareal_data_coarse_v(planeDataConfig),
		_parareal_data_coarse_u_prev(planeDataConfig), _parareal_data_coarse_v_prev(planeDataConfig),
		_parareal_data_output_u(planeDataConfig), _parareal_data_output_v(planeDataConfig),
		_parareal_data_output_u_prev(planeDataConfig), _parareal_data_output_v_prev(planeDataConfig),
		_parareal_data_error_u(planeDataConfig), _parareal_data_error_v(planeDataConfig)
#endif
	{
		// Calls initialization of the run (e.g. sets u, v)
		reset();

#if SWEET_PARAREAL
		parareal_setup();

#endif
	}

	virtual ~SimulationInstance()
	{
	}


	void reset()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "reset()" << std::endl;

		if (simVars.setup.benchmark_scenario_id <0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			BurgersValidationBenchmarks::printScenarioInformation();
			FatalError("Benchmark scenario not selected");
		}
		// Initialize diagnostics
		last_timestep_nr_update_diagnostics = -1;

		benchmark_diff_u = 0;
		benchmark_diff_v = 0;

		//TODO: is there a reason, why this is not called in swe_rexi
		simVars.reset();

		// set to some values for first touch NUMA policy (HPC stuff)

		tmp.physical_set_all(0);
		//Setup prog vars
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (param_use_staggering && simVars.disc.use_spectral_basis_diffs)
		{
			std::cerr << "Staggering and spectral basis not supported!" << std::endl;
			exit(1);
		}

		if (param_use_staggering && param_compute_error)
			std::cerr << "Warning: Staggered data will be interpolated to/from A-grid for exact linear solution" << std::endl;

		if (param_use_staggering)
		{
			/*
			 *              ^
			 *              |
			 *       ______v0,1_____
			 *       |             |
			 *       |			   |
			 *       |   (0.5,0.5) |
			 *  u0,0 |->  H/P0,0   |u1,0 ->
			 *(0,0.5)|	           |
			 *       |      ^      |
			 *   q0,0|______|______|
			 * (0,0)      v0,0
			 *           (0.5,0)
			 *
			 *
			 * These staggering should be used when interpolating from a staggered variable
			 * If interpolating from A grid to C staggered, use negative of displacements.
			 *
			 */
			stag_displacement[0] = 0.0;  // u_dx
			stag_displacement[1] = -0.5; // u_dy
			stag_displacement[2] = -0.5; // v_dx
			stag_displacement[3] = 0.0;  // v_dy
			stag_h[0] = -0.5;
			stag_h[1] = -0.5;
			stag_u[0] = -0.0;
			stag_u[1] = -0.5;
			stag_v[0] = -0.5;
			stag_v[1] = -0.0;
		}

		// Setup sampler for future interpolations
		sampler2D.setup(simVars.sim.domain_size, planeDataConfig);

		// Setup semi-Lagrangian
		semiLagrangian.setup(simVars.sim.domain_size, planeDataConfig);

		// Setup general (x,y) grid with position points

		PlaneData tmp_x(planeDataConfig);
		tmp_x.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = (double)i*(double)simVars.sim.domain_size[0]/(double)simVars.disc.res_physical[0];
			}
		);

		PlaneData tmp_y(planeDataConfig);
		tmp_y.physical_update_lambda_array_indices(
			[&](int i, int j, double &io_data)
			{
				io_data = ((double)j)*(double)simVars.sim.domain_size[1]/(double)simVars.disc.res_physical[1];
			}
		);

		pos_x = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_x);
		pos_y = Convert_PlaneData_To_ScalarDataArray::physical_convert(tmp_y);

		// Initialize arrival points with h position
		posx_a = pos_x+0.5*simVars.disc.cell_size[0];
		posy_a = pos_y+0.5*simVars.disc.cell_size[1];

		if (param_use_staggering)
		{
			prog_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);
			prog_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
				}
			);
		}
		else
		{
			prog_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);

			prog_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
				}
			);
		}

		//Initialize t-dt time step with initial condition
		prog_u_prev = prog_u;
		prog_v_prev = prog_v;

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
		{
			prog_u.file_physical_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		if (simVars.setup.input_data_filenames.size() > 1)
		{
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);
			if (param_use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

	// Print output info (if gui is disabled, this is done in main
		if (simVars.misc.gui_enabled)
			timestep_output();
	}


	// Calculate the model diagnostics
	void update_diagnostics()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "update_diagnostics()" << std::endl;

		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;


		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res_physical[0]*(double)simVars.disc.res_physical[1]);

		// energy
		simVars.diag.total_energy =
			0.5*((
					prog_u*prog_u +
					prog_v*prog_v
				).reduce_sum_quad()) * normalization;
	}


	/**
	 * Compute derivative for time stepping and store it to
	 * u_t and v_t
	 */
	void p_run_euler_timestep_update(
			const PlaneData &i_tmp,	///< prognostic variables
			const PlaneData &i_u,	///< prognostic variables
			const PlaneData &i_v,	///< prognostic variables

			PlaneData &o_tmp_t,		///< time updates
			PlaneData &o_u_t,		///< time updates
			PlaneData &o_v_t,		///< time updates

			double &o_dt,				///< time step restriction
			double i_fixed_dt = 0,		///< if this value is not equal to 0, use this time step size instead of computing one
			double i_simulation_timestamp = -1
	)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "p_run_euler_timestep_update()" << std::endl;

		/*
		 * 2D Burgers equation [with source term]
		 * u_t + u*u_x + v*u_y = nu*(u_xx + u_yy) [+ f(t,x,y)]
		 * v_t + u*v_x + v*v_y = nu*(v_xx + v_yy) [+ g(t,x,y)]
		 */

		//TODO: staggering vs. non staggering

		PlaneData f(planeDataConfig);
		BurgersValidationBenchmarks::set_source(i_simulation_timestamp,simVars,param_use_staggering,f);
		o_tmp_t.physical_set_all(0);

		/*
		 * u and v updates
		 */

		if (param_semilagrangian)
		{
			o_u_t = simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
			// Delete this line if no source is used.
			o_u_t += f;
			o_v_t = simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));
		}
		else
		{
			o_u_t = -(i_u*op.diff_c_x(i_u) + i_v*op.diff_c_y(i_u));
			o_u_t += simVars.sim.viscosity*(op.diff2_c_x(i_u)+op.diff2_c_y(i_u));
			// Delete this line if no source is used.
			o_u_t += f;
			o_v_t = -(i_u*op.diff_c_x(i_v) + i_v*op.diff_c_y(i_v));
			o_v_t += simVars.sim.viscosity*(op.diff2_c_x(i_v)+op.diff2_c_y(i_v));
		}

		o_tmp_t.physical_set_all(0);

		/*
		 * TIME STEP SIZE
		 */
		if (i_fixed_dt > 0)
		{
			o_dt = i_fixed_dt;
		}
		else
		{
			std::cout << "Only fixed time step size supported" << std::endl;
			assert(false);
			exit(1);
		}

	}

	/**
	 * IMEX time stepping for the coarse timestepping method
	 *
	 * The IMEX RK schemes combine an implicit and explicit time discretization.
	 * The diffusive, stiff term gets treated implicitly and the convective, non-
	 * stiff term gets treated explicitly. This results in the following system,
	 * which is solved by this routine:
	 * (I-\nu\Delta) u(t+\tau) = u(t) - \tau (u(t)*\nabla(u(t))
	 * for u(t+\tau).
	 */
	bool run_timestep_imex(
			PlaneData &io_u,
			PlaneData &io_v,

			double& o_dt,			///< return time step size for the computed time step
			double i_timestep_size,	///< timestep size

			PlaneOperators &op,
			const SimulationVariables &i_simVars,

			double i_max_simulation_time = std::numeric_limits<double>::infinity()	///< limit the maximum simulation time
			)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "run_timestep_imex()" << std::endl;

		PlaneData u=io_u;
		PlaneData v=io_v;

		// Modify timestep to final time if necessary
		double& t = o_dt;
		if (simVars.timecontrol.current_simulation_time+i_timestep_size < i_max_simulation_time)
			t = i_timestep_size;
		else
			t = i_max_simulation_time-simVars.timecontrol.current_simulation_time;

		// Initialize and set timestep dependent source for manufactured solution
		PlaneData f(planeDataConfig);
      PlaneData ff(planeDataConfig);
#if 0
      if (param_semilagrangian)
      {
		   BurgersValidationBenchmarks::set_source(simVars.timecontrol.current_simulation_time+t,simVars,param_use_staggering,f);
      }else
#endif
      {
		   BurgersValidationBenchmarks::set_source(simVars.timecontrol.current_simulation_time,simVars,param_use_staggering,f);
		   BurgersValidationBenchmarks::set_source(simVars.timecontrol.current_simulation_time+0.5*t,simVars,param_use_staggering,ff);
      }
		f.request_data_spectral();
      ff.request_data_spectral();

		// Setting explicit right hand side and operator of the left hand side
		PlaneData rhs_u = u;
		PlaneData rhs_v = v;

		if (param_semilagrangian)
		{
			rhs_u += t*f;
		}else{
			rhs_u += - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + 0.5*t*f;
			rhs_v += - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));
		}

		if (simVars.disc.use_spectral_basis_diffs) //spectral
		{
			PlaneData lhs = u;
			if (param_semilagrangian)
			{
#if 0
				lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
            io_u = rhs_u.spectral_div_element_wise(lhs);
            io_v = rhs_v.spectral_div_element_wise(lhs);
#else
//            std::cout << "*****************Warning******************" << std::endl << "explicit RK instead of SL!!!" << std::endl;
            PlaneData u1 = u + t*simVars.sim.viscosity*(op.diff2_c_x(u)+op.diff2_c_y(u))
                           - 0.5*t*(u*op.diff_c_x(u)+v*op.diff_c_y(u)) + f*t;
            PlaneData v1 = v + t*simVars.sim.viscosity*(op.diff2_c_x(v)+op.diff2_c_y(v))
                           - 0.5*t*(u*op.diff_c_x(v)+v*op.diff_c_y(v));

            io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
                  - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
            io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
                  - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));
#endif
			}else{
				lhs = ((-t)*simVars.sim.viscosity*(op.diff2_c_x + op.diff2_c_y)).spectral_addScalarAll(1.0);
            PlaneData u1 = rhs_u.spectral_div_element_wise(lhs);
            PlaneData v1 = rhs_v.spectral_div_element_wise(lhs);

            io_u = u + t*simVars.sim.viscosity*(op.diff2_c_x(u1)+op.diff2_c_y(u1))
                  - t*(u1*op.diff_c_x(u1)+v1*op.diff_c_y(u1)) +ff*t;
            io_v = v + t*simVars.sim.viscosity*(op.diff2_c_x(v1)+op.diff2_c_y(v1))
                  - t*(u1*op.diff_c_x(v1)+v1*op.diff_c_y(v1));
			}

		} else { //Jacobi
			FatalError("NOT available");
		}

		return true;
	}

	void run_timestep()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "run_timestep()" << std::endl;

		double dt = 0.0;

		// Only fixed time stepping supported with the Burgers equation
		assert(simVars.sim.CFL < 0);
		simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

#if 0
		if (param_semilagrangian)
		{
			dt = -simVars.sim.CFL;
			//Padding for last time step
			if (simVars.timecontrol.current_simulation_time+dt > simVars.timecontrol.max_simulation_time)
				dt = simVars.timecontrol.max_simulation_time-simVars.timecontrol.current_simulation_time;

			ScalarDataArray posx_d_tmp(planeDataConfig->physical_array_data_number_of_elements);
			ScalarDataArray posy_d_tmp(planeDataConfig->physical_array_data_number_of_elements);

			//Calculate departure points
			semiLagrangian.semi_lag_departure_points_settls(
							prog_u_prev, prog_v_prev,
							prog_u, prog_v,
							posx_a, posy_a,
							dt,
							posx_d, posy_d,
							stag_displacement
					);

			// Save old velocities
			prog_u_prev = prog_u;
			prog_v_prev = prog_v;

			//Now interpolate to the the departure points
			//Departure points are set for physical space
			prog_u = sampler2D.bicubic_scalar(
					prog_u,
					posx_d,
					posy_d,
					stag_u[0],
					stag_u[1]
			);

			prog_v = sampler2D.bicubic_scalar(
					prog_v,
					posx_d,
					posy_d,
					stag_v[0],
					stag_v[1]
			);

			if (simVars.disc.timestepping_method == 1)
			{
				// setup dummy data
				tmp.physical_set_all(0);

				// run standard Runge Kutta
				timestepping.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						tmp, prog_u, prog_v, ///< tmp is used to make use of the swe version of run_timestep
						dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
			}
			else if (simVars.disc.timestepping_method == 4)
			{
				// run IMEX Runge Kutta
				run_timestep_imex(
						prog_u, prog_v,
						dt,
						simVars.timecontrol.current_timestep_size,
						op,
						simVars,
						simVars.timecontrol.max_simulation_time
				);
			}
			else
				FatalError("Chosen time stepping method not available!");
		}
		else
#endif
		{
			if (simVars.disc.timestepping_method == 1)
			{
				// setup dummy data
				tmp.physical_set_all(0);

				// run standard Runge Kutta
				timestepping.run_timestep(
						this,
						&SimulationInstance::p_run_euler_timestep_update,	///< pointer to function to compute euler time step updates
						tmp, prog_u, prog_v, ///< tmp is used to make use of the swe version of run_timestep
						dt,
						simVars.timecontrol.current_timestep_size,
						simVars.disc.timestepping_order,
						simVars.timecontrol.current_simulation_time,
						simVars.timecontrol.max_simulation_time
					);
			}
			else if (simVars.disc.timestepping_method == 4)
			{
				// run IMEX Runge Kutta
				run_timestep_imex(
						prog_u, prog_v,
						dt,
						simVars.timecontrol.current_timestep_size,
						op,
						simVars,
						simVars.timecontrol.max_simulation_time
				);
			}
			else
				FatalError("Chosen time stepping method not available!");
		}

		//dt = simVars.timecontrol.current_timestep_size;

		// provide information to parameters
		simVars.timecontrol.current_timestep_size = dt;
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_output();
#endif
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const PlaneData &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "write_file()" << std::endl;

		char buffer[1024];

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		i_planeData.file_physical_saveData_ascii(buffer);

		return buffer;
	}


public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "timestep_output()" << std::endl;

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();

			// Print header
			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				o_ostream << "TIME\t\t\tTOT_ENERGY";
				if (param_compute_error){
					o_ostream << "\tMAX_ABS_U\tMAX_RMS_U\tMAX_U";
				}

				o_ostream << std::endl;

			}

			// Print timestep data to given output stream
			o_ostream << std::setprecision(8) << std::fixed << simVars.timecontrol.current_simulation_time << "\t" << simVars.diag.total_energy;

			if (param_compute_error){
				compute_errors(prog_u, prog_v);
				o_ostream << std::setprecision(8) << "\t" << benchmark_analytical_error_maxabs_u << "\t" << benchmark_analytical_error_rms_u << "\t" << prog_u.reduce_max();
			}

			o_ostream << std::endl;
		}

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			if (simVars.timecontrol.current_simulation_time != simVars.timecontrol.max_simulation_time)
				return false;

		// Dump data in csv, if requested
		if (simVars.misc.output_file_name_prefix.size() > 0)
		{
			write_file(prog_u, "u");
			write_file(prog_v, "v");
			if (param_compute_error)
				write_file(benchmark_analytical_error, "error");
		}

		if (simVars.misc.output_each_sim_seconds > 0)
			while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
				simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

		return true;
	}


	/*
	 * Compare manufactured solution with numerical solution
	 */
public:
	void compute_errors(
         const PlaneData &i_planeData_u,
         const PlaneData &i_planeData_v
	)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "compute_errors()" << std::endl;

		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData u = i_planeData_u;
		PlaneData v = i_planeData_v;

		// Only possible for manufactured solutions
		if (simVars.setup.benchmark_scenario_id < 51 && simVars.setup.benchmark_scenario_id > 59)
			return;

		//Analytical solution at current time on original grid (stag or not)
		PlaneData ts_u(planeDataConfig);
		PlaneData ts_v(planeDataConfig);

		if (param_use_staggering)
		{
			ts_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);

			ts_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];
					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
				}
			);
		}
		else
		{
			ts_u.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);

			ts_v.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.res_physical[0])*simVars.sim.domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.res_physical[1])*simVars.sim.domain_size[1];

					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
				}
			);
		}


		benchmark_analytical_error = ts_u-u;

		benchmark_analytical_error_rms_u = (ts_u-u).reduce_rms_quad();
		benchmark_analytical_error_rms_v = (ts_v-v).reduce_rms_quad();

		benchmark_analytical_error_maxabs_u = (ts_u-u).reduce_maxAbs();
		benchmark_analytical_error_maxabs_v = (ts_v-v).reduce_maxAbs();
	}


	bool should_quit()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "should_quit()" << std::endl;

		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		if (simVars.timecontrol.max_simulation_time != -1 && simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time)
			return true;

		return false;
	}



	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_post_frame_processing()" << std::endl;

		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}



	struct VisStuff
	{
		const PlaneData* data;
		const char *description;
	};

	/**
	 * Arrays for online visualisation and their textual description
	 */
	VisStuff vis_arrays[NUM_OF_UNKNOWNS] =
	{
			{&prog_u,	"u"},
			{&prog_v,	"v"}
	};

	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
	)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_get_vis_data_array()" << std::endl;

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
		*o_dataArray = vis_arrays[id].data;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
	}


	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_get_status_string()" << std::endl;

		// first, update diagnostic values if required
		update_diagnostics();

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		static char title_string[1024];
		sprintf(title_string, "Time: %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %.14s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				vis_arrays[id].description,
				simVars.diag.total_mass, simVars.diag.total_energy, simVars.diag.total_potential_enstrophy);
		return title_string;
	}



	void vis_pause()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_pause()" << std::endl;

		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_keypress()" << std::endl;

		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			if (simVars.misc.vis_id >= NUM_OF_UNKNOWNS)
				simVars.misc.vis_id = 0;
			break;

		case 'V':
			simVars.misc.vis_id--;
			if (simVars.misc.vis_id < 0)
				simVars.misc.vis_id = 1;
			break;
		}
	}


	bool instability_detected()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "instability_detected()" << std::endl;

		// Necessary to circumvent FFTW transformations on prog_u and prog_v, which would lead to errors
		PlaneData u = prog_u;
		PlaneData v = prog_v;

		return !(u.reduce_boolean_all_finite() && v.reduce_boolean_all_finite());
	}


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	PlaneData _parareal_data_start_u, _parareal_data_start_v,
		_parareal_data_start_u_prev, _parareal_data_start_v_prev;
	Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2> parareal_data_start;

	PlaneData _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_PlaneData<NUM_OF_UNKNOWNS> parareal_data_fine;

	PlaneData _parareal_data_coarse_u, _parareal_data_coarse_v,
		_parareal_data_coarse_u_prev, _parareal_data_coarse_v_prev;
	Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2> parareal_data_coarse;

	PlaneData _parareal_data_output_u, _parareal_data_output_v,
		_parareal_data_output_u_prev, _parareal_data_output_v_prev;
	Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2> parareal_data_output;

	PlaneData _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_PlaneData<NUM_OF_UNKNOWNS> parareal_data_error;

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "parareal_setup()" << std::endl;

		{
			PlaneData* data_array[NUM_OF_UNKNOWNS*2] = {&_parareal_data_start_u, &_parareal_data_start_v,
					&_parareal_data_start_u_prev, &_parareal_data_start_v_prev};
			parareal_data_start.setup(data_array);
		}

		{
			PlaneData* data_array[NUM_OF_UNKNOWNS] = {&_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			PlaneData* data_array[NUM_OF_UNKNOWNS*2] = {&_parareal_data_coarse_u, &_parareal_data_coarse_v,
					&_parareal_data_coarse_u_prev, &_parareal_data_coarse_v_prev};
			parareal_data_coarse.setup(data_array);
		}

		{
			PlaneData* data_array[NUM_OF_UNKNOWNS*2] = {&_parareal_data_output_u, &_parareal_data_output_v,
					&_parareal_data_output_u_prev, &_parareal_data_output_v_prev};
			parareal_data_output.setup(data_array);
		}

		{
			PlaneData* data_array[NUM_OF_UNKNOWNS] = {&_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

		output_data_valid = false;
	}



	/**
	 * Set the start and end of the coarse time step
	 */
	void sim_set_timeframe(
			double i_timeframe_start,	///< start timestamp of coarse time step
			double i_timeframe_end		///< end time stamp of coarse time step
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "Timeframe: [" << i_timeframe_start << ", " << i_timeframe_end << "]" << std::endl;

		timeframe_start = i_timeframe_start;
		timeframe_end = i_timeframe_end;
	}



	/**
	 * Set the initial data at i_timeframe_start
	 */
	void sim_setup_initial_data(
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_setup_initial_data()" << std::endl;

		reset();


		*parareal_data_start.data_arrays[0] = prog_u;
		*parareal_data_start.data_arrays[1] = prog_v;
		*parareal_data_start.data_arrays[2] = prog_u_prev;
		*parareal_data_start.data_arrays[3] = prog_v_prev;

	}

	/**
	 * Set simulation data to data given in i_sim_data.
	 * This can be data which is computed by another simulation.
	 * Y^S := i_sim_data
	 */
	void sim_set_data(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data()" << std::endl;

		// copy to buffers
		parareal_data_start = i_pararealData;

		// cast to pararealPlaneData stuff
	}

	/**
	 * Set the MPI communicator to use for simulation purpose
	 * (TODO: not yet implemented since our parallelization-in-space
	 * is done only via OpenMP)
	 */
	void sim_set_mpi_comm(
			int i_mpi_comm
	)
	{
		// NOTHING TO DO HERE
	}

	/**
	 * compute solution on time slice with fine timestep:
	 * Y^F := F(Y^S)
	 */
	void run_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_fine()" << std::endl;

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		bool was_sl = false;
		if (param_semilagrangian)
		{
			param_semilagrangian = false;
			was_sl = true;
		}

		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		if (was_sl)
			param_semilagrangian = true;

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_u;
		*parareal_data_fine.data_arrays[1] = prog_v;
	}


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_Data& get_reference_to_data_timestep_fine()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_fine()" << std::endl;

		return parareal_data_fine;
	}

#endif

#if SWEET_PARAREAL

	/**
	 * compute solution with coarse timestepping:
	 * Y^C := G(Y^S)
	 */
	void run_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "run_timestep_coarse()" << std::endl;

		prog_u = *parareal_data_start.data_arrays[0];
		prog_v = *parareal_data_start.data_arrays[1];
		prog_u_prev = *parareal_data_start.data_arrays[2];
		prog_v_prev = *parareal_data_start.data_arrays[3];

		//Preserve parareal_data_start for next timestep to be prog_u_prev
		*parareal_data_output.data_arrays[2] = prog_u;
		*parareal_data_output.data_arrays[3] = prog_v;
		*parareal_data_coarse.data_arrays[2] = prog_u;
		*parareal_data_coarse.data_arrays[3] = prog_v;

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// set Runge-Kutta scheme to the chosen one for coarse time stepping
		int tmpScheme = simVars.disc.timestepping_method;
		int tmpOrder = simVars.disc.timestepping_order;
		simVars.disc.timestepping_method = simVars.disc.timestepping_method2;
		simVars.disc.timestepping_order = simVars.disc.timestepping_order2;
		// save the fine delta t to restore it later
		double tmpCFL = simVars.sim.CFL;
		simVars.sim.CFL=timeframe_start-timeframe_end;

		// make multiple time steps in the coarse solver possible
		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// restore fine delta t and time stepping scheme
		simVars.sim.CFL = tmpCFL;
		simVars.disc.timestepping_method = tmpScheme;
		simVars.disc.timestepping_order = tmpOrder;
		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_u;
		*parareal_data_coarse.data_arrays[1] = prog_v;
	}



	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_Data& get_reference_to_data_timestep_coarse()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_data_timestep_coarse()" << std::endl;

		return parareal_data_coarse;
	}



	/**
	 * Compute the error between the fine and coarse timestepping:
	 * Y^E := Y^F - Y^C
	 */
	void compute_difference()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_difference()" << std::endl;

		for (int k = 0; k < 2; k++)
			*parareal_data_error.data_arrays[k] = *parareal_data_fine.data_arrays[k] - *parareal_data_coarse.data_arrays[k];
	}



	/**
	 * Compute the data to be forwarded to the next time step
	 * Y^O := Y^C + Y^E
	 *
	 * Return: Error indicator based on the computed error norm between the
	 * old values and new values
	 */
	double compute_output_data(
			bool i_compute_convergence_test
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "compute_output_data()" << std::endl;

		double convergence = -1;

		if (!i_compute_convergence_test || !output_data_valid)
		{
			for (int k = 0; k < NUM_OF_UNKNOWNS; k++)
				*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			output_data_valid = true;
			return convergence;
		}



		for (int k = 0; k < NUM_OF_UNKNOWNS; k++)
		{
			tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_u = *parareal_data_output.data_arrays[0];
		prog_v = *parareal_data_output.data_arrays[1];

		if (param_compute_error){
			compute_errors(prog_u, prog_v);
			std::cout << "maxabs error compared to analytical solution: " << benchmark_analytical_error_maxabs_u << std::endl;
		}

		output_data_valid = true;
		return convergence;
	}



	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_Data& get_reference_to_output_data()
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "get_output_data()" << std::endl;

		return parareal_data_output;
	}


	void output_data_file(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "output_data_file()" << std::endl;

		Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2>& data = (Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2>&)i_data;

		std::ostringstream ss;
		ss << simVars.misc.output_file_name_prefix << "_iter" << iteration_id << "_slice" << time_slice_id << ".csv";

		std::string filename = ss.str();

		data.data_arrays[0]->file_physical_saveData_ascii(filename.c_str());
		//data.data_arrays[0]->file_saveSpectralData_ascii(filename.c_str());

		std::ostringstream ss2;
		ss2 << simVars.misc.output_file_name_prefix << "_iter" << iteration_id << "_slice" << time_slice_id << ".err";

		std::string filename2 = ss2.str();

		compute_errors(*data.data_arrays[0], *data.data_arrays[1]);

		benchmark_analytical_error.file_physical_saveData_ascii(filename2.c_str());

		//data.data_arrays[0]->file_saveData_vtk(filename.c_str(), filename.c_str());
	}


	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "output_data_console()" << std::endl;

		update_diagnostics();
		// Print timestep data to console
		std::cout << std::setprecision(8) << "Total energy: " << simVars.diag.total_energy << std::endl;
	}

#endif
};




int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

	MemBlockAlloc::setup();

	// program specific input parameter names
	const char *bogus_var_names[] = {
			"compute-error",
			"staggering",
			"semi-lagrangian",
			nullptr
	};

	// default values for program specific parameter
	simVars.bogus.var[0] = param_compute_error; // compute-error
	simVars.bogus.var[1] = param_use_staggering; // staggering
	simVars.bogus.var[2] = param_semilagrangian; // semi-lagrangian

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--compute-error [0/1]		Compute the errors, default=0" << std::endl;
		std::cout << "	--staggering [0/1]		Use staggered grid, default=0" << std::endl;
		std::cout << "	--semi-lagrangian [0/1]		Use semi-lagrangian formulation, default=0" << std::endl;
		std::cout << std::endl;

#if SWEET_PARAREAL
		simVars.parareal.setup_printOptions();
#endif
		return -1;
	}

	//Burgers parameters
	param_compute_error = simVars.bogus.var[0];
	param_use_staggering = simVars.bogus.var[1];
	param_semilagrangian = simVars.bogus.var[2];

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	std::cout << planeDataConfigInstance.getUniqueIDString() << std::endl;


	std::ostringstream buf;
	buf << std::setprecision(14);

	{

#if SWEET_PARAREAL
		if (simVars.parareal.enabled)
		{
			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller_Serial<SimulationInstance> parareal_Controller_Serial;

			// setup controller. This initializes several simulation instances
			parareal_Controller_Serial.setup(&simVars.parareal);

			// execute the simulation
			parareal_Controller_Serial.run();
		}
		else
#endif

#if SWEET_GUI
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationBurgers = new SimulationInstance;
			VisSweet<SimulationInstance> visSweet(simulationBurgers);
			delete simulationBurgers;
		}
		else
#endif
		{
			SimulationInstance *simulationBurgers = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			simulationBurgers->reset();

			//Time counter
			Stopwatch time;

			//Diagnostic measures at initial stage
			double diagnostics_energy_start, diagnostics_mass_start, diagnostics_potential_entrophy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
				simulationBurgers->update_diagnostics();
				diagnostics_energy_start = simVars.diag.total_energy;
				diagnostics_mass_start = simVars.diag.total_mass;
				diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;
			}

			//Start counting time
			time.reset();

			//Main time loop
			while(true)
			{
				//Output data
				if (simulationBurgers->timestep_output(buf))
				{
					// string output data
					std::string output = buf.str();
					buf.str("");

					// This is an output printed on screen or buffered to files if > used
					std::cout << output;
				}

				//Stop simulation if requested
				if (simulationBurgers->should_quit())
					break;

				//Main call for timestep run
				simulationBurgers->run_timestep();

				//Instability
				if (simulationBurgers->instability_detected())
				{
					std::cout << "INSTABILITY DETECTED" << std::endl;
					break;
				}
			}

			//Stop counting time
			time.stop();

			double seconds = time();

			//End of output results
			std::cout << "Simulation time (seconds): " << seconds << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << seconds/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;


			if (simVars.misc.verbosity > 0)
			{
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs(simVars.diag.total_energy-diagnostics_energy_start) << std::endl;
				std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs(simVars.diag.total_mass-diagnostics_mass_start) << std::endl;
				std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs(simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start) << std::endl;

				if (simVars.setup.benchmark_scenario_id == 2 || simVars.setup.benchmark_scenario_id == 3 || simVars.setup.benchmark_scenario_id == 4)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationBurgers->benchmark_diff_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationBurgers->benchmark_diff_v << std::endl;
				}
			}

			delete simulationBurgers;
		}
	}

	return 0;
}
