/*
 * Burgers equation
 *
 */

#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataSemiLagrangian.hpp>
#include <sweet/plane/Convert_ScalarDataArray_to_PlaneData.hpp>
#include <sweet/FatalError.hpp>

#include "burgers/Burgers_Plane_TimeSteppers.hpp"

#include <sweet/Stopwatch.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal.hpp>
#endif


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

// Input parameters (cmd line)

// General parameters
SimulationVariables simVars;

const int NUM_OF_UNKNOWNS=2;

//specific parameters
bool param_compute_error = 1;


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

	// implementation of different time steppers
	Burgers_Plane_TimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// implementation of different time steppers
	Burgers_Plane_TimeSteppers timeSteppersCoarse;
#endif

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
	PlaneDataTimesteppingRK timestepping_rk;


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

		// Variable dimensions (mem. allocation)
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		prog_u_prev(planeDataConfig),
		prog_v_prev(planeDataConfig),

		benchmark_analytical_error(planeDataConfig),

		// Initialise operators
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
		if (simVars.parareal.enabled)
			parareal_setup();

#endif
	}


	virtual ~SimulationInstance()
	{
	}


	/*
	 * Reset all variables
	 */
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
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_u.spectral_set_all(0,0);
		prog_v.spectral_set_all(0,0);
#endif
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (simVars.disc.use_staggering && simVars.disc.use_spectral_basis_diffs)
		{
			std::cerr << "Staggering and spectral basis not supported!" << std::endl;
			exit(1);
		}

		if (simVars.disc.use_staggering && param_compute_error)
			std::cerr << "Warning: Staggered data will be interpolated to/from A-grid for exact linear solution" << std::endl;


		// Set initial condigions given from BurgersValidationBenchmarks
		if (simVars.disc.use_staggering)
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
			if (simVars.disc.use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		if (simVars.setup.input_data_filenames.size() > 1)
		{
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);
			if (simVars.disc.use_staggering && param_compute_error)
			{
				std::cerr << "Warning: Computing analytical solution with staggered grid based on loaded data not supported!" << std::endl;
				exit(1);
			}
		}

		timeSteppers.setup(
			simVars.disc.timestepping_method,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			op,
			simVars
		);

		timestep_output();
	}


	/*
	 * Calculate the model diagnostics
	 */
	void update_diagnostics()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "update_diagnostics()" << std::endl;

		// Assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;

		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res_physical[0]*(double)simVars.disc.res_physical[1]);

		// Energy
		simVars.diag.total_energy =
			0.5*((
					prog_u*prog_u +
					prog_v*prog_v
				).reduce_sum_quad()) * normalization;
	}


	/*
	 * Routine to do one time step of chosen scheme and order
	 */
	void run_timestep()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "run_timestep()" << std::endl;

		double dt = 0.0;

		// Only fixed time stepping supported with the Burgers equation
		assert(simVars.sim.CFL < 0);
		simVars.timecontrol.current_timestep_size = -simVars.sim.CFL;

		timeSteppers.master->run_timestep(
			prog_u, prog_v,
			prog_u_prev, prog_v_prev,
			dt,
			simVars.timecontrol.current_timestep_size,
			simVars.timecontrol.current_simulation_time,
			simVars.timecontrol.max_simulation_time
		);

#if 0
		if (simVars.disc.timestepping_method == "ln_erk") //Explicit RK
		{
			// Basic explicit Runge-Kutta
			if (simVars.disc.timestepping_order != 21)
			{
				// setup dummy data
				tmp.physical_set_all(0);

				// run standard Runge Kutta
				timestepping_rk.run_timestep(
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
			// Explicit Runge-Kutta with order 1 in diffusion and order 2 in advection
			else
			{
				burgers_plane.run_timestep_explicit_ts(
						prog_u,
						prog_v,

						dt,
						simVars.timecontrol.current_timestep_size,

						op,
						simVars,
						planeDataConfig
				);
			}
		}
		else if (simVars.disc.timestepping_method == "ln_imex") //IMEX
		{
			burgers_plane.run_timestep_imex(
					prog_u, prog_v,
					dt,
					simVars.timecontrol.current_timestep_size,
					op,
					simVars,
					planeDataConfig,
					simVars.timecontrol.max_simulation_time
			);
		}
		else if (simVars.disc.timestepping_method == "l_irk_n_sl") //SL IMPLICIT
		{
			burgers_plane.run_timestep_sl(
				prog_u, prog_v,
				prog_u_prev, prog_v_prev,
				posx_a, posy_a,
				dt,
				simVars.timecontrol.current_timestep_size,
				simVars,
				planeDataConfig,
				op,
				sampler2D,
				semiLagrangian,
				staggering
			);
		}
		else
		{
			FatalError("Chosen time stepping method not available!");
		}

		//dt = simVars.timecontrol.current_timestep_size;
#endif
		// Provide information to parameters
		simVars.timecontrol.current_timestep_size = dt;
		simVars.timecontrol.current_simulation_time += dt;
		simVars.timecontrol.current_timestep_nr++;

		timestep_output();
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


	/*
	 * Create stream with output data
	 */
public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "timestep_output()" << std::endl;

		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

		if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.current_simulation_time)
			return false;

		if (param_compute_error)
			compute_errors(prog_u, prog_v);

		// Dump data in csv, if requested
		if (simVars.misc.output_file_name_prefix.size() > 0)
		{
			write_file(prog_u, "prog_u");
			write_file(prog_v, "prog_v");
			if (param_compute_error)
				write_file(benchmark_analytical_error, "error");
		}

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

			if (param_compute_error)
				o_ostream << std::setprecision(8) << "\t" << benchmark_analytical_error_maxabs_u << "\t" << benchmark_analytical_error_rms_u << "\t" << prog_u.reduce_max();

			o_ostream << std::endl;
		}

		if (simVars.misc.output_each_sim_seconds > 0)
		{
			if (simVars.misc.output_next_sim_seconds == simVars.timecontrol.max_simulation_time)
			{
				simVars.misc.output_next_sim_seconds = std::numeric_limits<double>::infinity();
			}
			else
			{
				while (simVars.misc.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

				if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.max_simulation_time)
					simVars.misc.output_next_sim_seconds = simVars.timecontrol.max_simulation_time;
			}
		}

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

		// Analytical solution at current time on original grid (stag or not)
		PlaneData ts_u(planeDataConfig);
		PlaneData ts_v(planeDataConfig);

		if (simVars.disc.use_staggering)
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


	/*
	 * Check for final time step
	 */
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

#if SWEET_GUI

	/**
	 * Postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_post_frame_processing()" << std::endl;

		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations; i++)
				run_timestep();
	}


	/*
	 * Struct for visualization
	 */
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


	/*
	 * Getter for visualization data array
	 */
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
	 * Return status string for window title
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


/*
 * Pause visualization
 */
	void vis_pause()
	{
		if (simVars.misc.verbosity > 2)
			std::cout << "vis_pause()" << std::endl;

		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}


/*
 * Handle keypress events in visualization
 */
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
#endif


	/*
	 * Detect instabilities in the calculation
	 */
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


	/*
	 * Setup Parareal variables
	 */
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

		timeSteppers.setup(
				simVars.disc.timestepping_method,
				simVars.disc.timestepping_order,
				simVars.disc.timestepping_order2,
				op,
				simVars
			);

		timeSteppersCoarse.setup(
				simVars.parareal.coarse_timestepping_method,
				simVars.parareal.coarse_timestepping_order,
				simVars.parareal.coarse_timestepping_order2,
				op,
				simVars
			);

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

		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			this->run_timestep(simVars.disc.timestepping_method, simVars.disc.timestepping_order);
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

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

		// Preserve parareal_data_start for next timestep to be prog_u_prev
		*parareal_data_output.data_arrays[2] = prog_u;
		*parareal_data_output.data_arrays[3] = prog_v;
		*parareal_data_coarse.data_arrays[2] = prog_u;
		*parareal_data_coarse.data_arrays[3] = prog_v;

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		simVars.timecontrol.current_timestep_size=timeframe_end-timeframe_start;
		double dt = 0.0;

		// make multiple time steps in the coarse solver possible
		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			//run_timestep(simVars.parareal.coarse_timestepping_method, simVars.parareal.coarse_timestepping_order);
			timeSteppersCoarse.master->run_timestep(
				prog_u, prog_v,
				prog_u_prev, prog_v_prev,
				dt,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.max_simulation_time
			);
			// Provide information to parameters
			simVars.timecontrol.current_timestep_size = dt;
			simVars.timecontrol.current_simulation_time += dt;
			simVars.timecontrol.current_timestep_nr++;

			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

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

		PlaneData tmp(planeDataConfig);
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		tmp.spectral_set_all(0,0);
#endif
		tmp.physical_set_all(0.0);

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


	/*
	 * Write output of Parareal data to file
	 */
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


	/*
	 * Write output of Parareal data to console
	 */
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
			nullptr
	};

	// default values for program specific parameter
	simVars.bogus.var[0] = param_compute_error; // compute-error

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--compute-error [0/1]		Compute errors (if available, default: 1)" << std::endl;
		std::cout << std::endl;

/*
 * #if SWEET_PARAREAL
 * 		simVars.parareal.setup_printOptions();
 * #endif
 */
		return -1;
	}

	//Burgers parameters
	param_compute_error = simVars.bogus.var[0];

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	// Print header
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Solving viscous Burgers' equation" << std::endl;
	std::cout << "-----------------------------" << std::endl;
	simVars.outputConfig();
	std::cout << "LOCAL PARAMETERS" << std::endl;
	std::cout << " + param_compute_error: " << param_compute_error << std::endl;
	std::cout << std::endl;

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

			// Setting initial conditions and workspace - in case there is no GUI
			simulationBurgers->reset();

			// Time counter
			Stopwatch time;

			// Diagnostic measures at initial stage
			double diagnostics_energy_start;

			// Initialize diagnostics
			if (simVars.misc.verbosity > 0)
			{
				simulationBurgers->update_diagnostics();
				diagnostics_energy_start = simVars.diag.total_energy;
			}

			// Start counting time
			time.reset();

			// Main time loop
			while(true)
			{
				/*
				//Output data
				if (simulationBurgers->timestep_output(buf))
				{
					// string output data
					std::string output = buf.str();
					buf.str("");

					// This is an output printed on screen or buffered to files if > used
					std::cout << output;
				}
				*/

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

				if (param_compute_error)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationBurgers->benchmark_analytical_error_maxabs_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationBurgers->benchmark_analytical_error_maxabs_v << std::endl;
				}
			}

			delete simulationBurgers;
		}
	}

	return 0;
}
