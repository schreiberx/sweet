/*
 * Burgers equation
 *
 */

#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/Convert_ScalarDataArray_to_PlaneData.hpp>



#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include "burgers/Burgers_Plane_TimeSteppers.hpp"



// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal.hpp>
#endif


// General parameters
SimulationVariables simVars;

const int NUM_OF_UNKNOWNS=2;


class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	// u: velocity in x-direction
	// v: velocity in y-direction
	PlaneData prog_u, prog_v;

	// Prognostic variables at time step t-dt
	PlaneData prog_u_prev, prog_v_prev;

#if SWEET_GUI
	// visualization variable
	PlaneData vis;
#endif

	// Initial values for comparison with analytical solution
	PlaneData t0_prog_u, t0_prog_v;

	// implementation of different time steppers
	Burgers_Plane_TimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// implementation of different time steppers
	Burgers_Plane_TimeSteppers timeSteppersCoarse;
#endif

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

   // Variable to select correct analytic solver
   int analytic_solution = 2;

	class BenchmarkErrors
	{
	public:
		// Max difference to initial conditions
		double benchmark_diff_u;
		double benchmark_diff_v;

		// Error measures L2 norm
		double benchmark_analytical_error_rms_u;
		double benchmark_analytical_error_rms_v;

		// Error measures max norm
		double benchmark_analytical_error_maxabs_u;
		double benchmark_analytical_error_maxabs_v;
	};

	BenchmarkErrors benchmark;

	// Finite difference operators
	PlaneOperators op;


	// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;

public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		prog_u_prev(planeDataConfig),
		prog_v_prev(planeDataConfig),

#if SWEET_GUI
		vis(planeDataConfig),
#endif

		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

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



	void reset()
	{
		if (simVars.setup.benchmark_scenario_id <0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			BurgersValidationBenchmarks::printScenarioInformation();
			FatalError("Benchmark scenario not selected");
		}

		// Initialize diagnostics
		last_timestep_nr_update_diagnostics = -1;

		benchmark.benchmark_diff_u = 0;
		benchmark.benchmark_diff_v = 0;
		benchmark.benchmark_analytical_error_rms_u = 0;
		benchmark.benchmark_analytical_error_rms_v = 0;
		benchmark.benchmark_analytical_error_maxabs_u = 0;
		benchmark.benchmark_analytical_error_maxabs_v = 0;

		//TODO: is there a reason, why this is not called in swe_rexi
		simVars.reset();

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_u.spectral_set_all(0,0);
		prog_v.spectral_set_all(0,0);
		prog_u_prev.spectral_set_all(0,0);
		prog_v_prev.spectral_set_all(0,0);
#endif
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);
		prog_u_prev.physical_set_all(0);
		prog_v_prev.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (simVars.disc.use_staggering && simVars.disc.use_spectral_basis_diffs)
			FatalError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.use_staggering ||  !simVars.disc.use_spectral_basis_diffs)
			FatalError("Finite differences and spectral dealiasing should not be used together! Please compile without dealiasing.");
#endif

		timeSteppers.setup(
			simVars.disc.timestepping_method,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			op,
			simVars
		);

		// Set initial conditions given from BurgersValidationBenchmarks
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

		// Setup look-up table of Manufactured Solution
		if (simVars.disc.timestepping_method.find("_mms")!=std::string::npos)
		{
#if SWEET_PARAREAL
			if (!simVars.parareal.enabled)
#endif
			timeSteppers.setup_look_up_table(
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.max_simulation_time,
				simVars.timecontrol.current_timestep_size
			);

			// Set before anyway
//			timeSteppers.return_initial(prog_u);
//#if SWEET_USE_PLANE_SPECTRAL_SPACE
//			prog_v.spectral_set_all(0,0);
//#endif
//			prog_v.physical_set_all(0);
		}

		// Initialize t-dt time step with initial condition
		prog_u_prev = prog_u;
		prog_v_prev = prog_v;
		// Save initial conditions for analytic solution
		t0_prog_u = prog_u;
		t0_prog_v = prog_v;

		// Load data, if requested
		if (simVars.setup.input_data_filenames.size() > 0)
			prog_u.file_physical_loadData(simVars.setup.input_data_filenames[0].c_str(), simVars.setup.input_data_binary);

		if (simVars.setup.input_data_filenames.size() > 1)
			prog_v.file_physical_loadData(simVars.setup.input_data_filenames[1].c_str(), simVars.setup.input_data_binary);


		update_diagnostics();
		diagnostics_energy_start = simVars.diag.total_energy;

#if SWEET_PARAREAL
		if (!simVars.parareal.enabled)
#endif
		timestep_output();

      if (simVars.misc.compute_errors)
      {
         bool foundl = (simVars.disc.timestepping_method.find("l_")==0) || (simVars.disc.timestepping_method.find("_l_")!=std::string::npos);
         bool foundn = (simVars.disc.timestepping_method.find("n_")==0) || (simVars.disc.timestepping_method.find("_n_")!=std::string::npos);
         bool foundnl = (simVars.disc.timestepping_method.find("ln_")==0) || (foundl && foundn);

         if (foundnl)
            analytic_solution = 1;
         else if (foundl)
            analytic_solution = 2;
         else
            FatalError("Computing errors for this timestepping-method is not possible");
      }
	}


	// Calculate the model diagnostics
	void update_diagnostics()
	{
#if !SWEET_PARAREAL
		// Assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;
#endif

		double normalization = (simVars.sim.domain_size[0]*simVars.sim.domain_size[1]) /
								((double)simVars.disc.res_physical[0]*(double)simVars.disc.res_physical[1]);

		// Reduce amount of possible FFTs to minimize numerical error
		PlaneData tmp_u = prog_u;
		PlaneData tmp_v = prog_v;

		// Energy
		simVars.diag.total_energy =
			0.5*((
					tmp_u*tmp_u +
					tmp_v*tmp_v
				).reduce_sum_quad()) * normalization;
	}


	/*
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		if (simVars.disc.timestepping_method == "ln_cole_hopf" || simVars.disc.timestepping_method == "l_direct")
		{
			prog_u = t0_prog_u;
			prog_v = t0_prog_v;
			timeSteppers.master->run_timestep(
				prog_u, prog_v,
				prog_u_prev, prog_v_prev,
				simVars.timecontrol.current_timestep_size+simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time
			);
		}
		else
		{
			timeSteppers.master->run_timestep(
				prog_u, prog_v,
				prog_u_prev, prog_v_prev,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);
		}

		// Advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			FatalError("Max simulation time exceeded!");

#if SWEET_PARAREAL
		if (!simVars.parareal.enabled)
#endif
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
		char buffer[1024];

		const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
//		i_planeData.file_physical_saveData_ascii(buffer);
		i_planeData.file_physical_saveData_ascii(buffer,'\n',16,1);

//		char tmp[128];
//		strcpy(tmp,i_name);
//		strcat(tmp,"_spec");
//		sprintf(buffer, filename_template, tmp, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
//		i_planeData.file_spectral_saveData_ascii(buffer,'\n',12,1);
//		i_planeData.file_spectral_abs_arg_saveData_ascii(buffer,'\n',12,1);

		return buffer;
	}



public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		// output each time step
		if (simVars.misc.output_each_sim_seconds < 0)
			return false;

//		if (simVars.misc.output_next_sim_seconds-simVars.misc.output_next_sim_seconds*(1e-12) > simVars.timecontrol.current_simulation_time)
//			return false;
		if (simVars.misc.output_next_sim_seconds-simVars.timecontrol.current_timestep_size*0.5 > simVars.timecontrol.current_simulation_time)
			return false;

		PlaneData tmp_u = prog_u;
		PlaneData tmp_v = prog_v;

		// Dump data in csv, if requested
		if (simVars.misc.output_file_name_prefix.size() > 0)
		{
			write_file(tmp_u, "prog_u");
			//write_file(tmp_v, "prog_v");

			tmp_u.request_data_spectral();

			char buffer[1024];
			const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
			sprintf(buffer,filename_template,"prog_u_amp_phase",simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);

			std::ofstream file(buffer, std::ios_base::trunc);
			file << std::setprecision(12);

			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				file << x << ", " << tmp_u.spectral_return_amplitude(0,x) << ", " << tmp_u.spectral_return_phase(0,x) << std::endl;
			}
			file.close();
		}

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();
			if (simVars.misc.compute_errors)
			{
				PlaneData tmp(planeDataConfig);
				tmp = compute_errors2(tmp_u,tmp_v);

				if (simVars.misc.output_file_name_prefix.size() > 0)
				{
					write_file(tmp, "analytical");

					tmp.request_data_spectral();

					char buffer[1024];
					const char* filename_template = simVars.misc.output_file_name_prefix.c_str();
					sprintf(buffer,filename_template,"analytical_amp_phase",simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);

					std::ofstream file(buffer, std::ios_base::trunc);
					file << std::setprecision(12);
					for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
					{
						file << x << ", " << tmp.spectral_return_amplitude(0,x) << ", " << tmp.spectral_return_phase(0,x) << std::endl;
					}
					file.close();
				}
			}

			std::stringstream header;
			std::stringstream rows;

			rows << std::setprecision(12);

			// Prefix
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "DATA";
			rows << "DATA";

			// Time
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "\tT";
			rows << "\t" << simVars.timecontrol.current_simulation_time;

			// Energy
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "\tTOTAL_ENERGY";
			rows << "\t" << simVars.diag.total_energy;

			if (simVars.misc.compute_errors)
			{
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tMAX_ABS_U\tMAX_RMS_U\tMAX_U";
				rows << "\t" << benchmark.benchmark_analytical_error_maxabs_u << "\t" << benchmark.benchmark_analytical_error_rms_u << "\t" << prog_u.reduce_max();
			}

			if (simVars.timecontrol.current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;
		}

		if (simVars.misc.output_each_sim_seconds > 0)
		{
			if (simVars.misc.output_next_sim_seconds == simVars.timecontrol.max_simulation_time)
			{
				simVars.misc.output_next_sim_seconds = std::numeric_limits<double>::infinity();
			}
			else
			{
//				while (simVars.misc.output_next_sim_seconds-simVars.misc.output_next_sim_seconds*(1e-12) <= simVars.timecontrol.current_simulation_time)
//					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;
				while (simVars.misc.output_next_sim_seconds-simVars.timecontrol.current_timestep_size*0.5 <= simVars.timecontrol.current_simulation_time)
					simVars.misc.output_next_sim_seconds += simVars.misc.output_each_sim_seconds;

				if (simVars.misc.output_next_sim_seconds > simVars.timecontrol.max_simulation_time)
					simVars.misc.output_next_sim_seconds = simVars.timecontrol.max_simulation_time;
			}
		}

		return true;
	}


public:
	void compute_errors(
         const PlaneData &i_planeData_u,
         const PlaneData &i_planeData_v
	)
	{
		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData u = i_planeData_u;
		PlaneData v = i_planeData_v;

		// Analytical solution at current time on original grid
		PlaneData ts_u = t0_prog_u;
		PlaneData ts_v = t0_prog_v;

		if (simVars.misc.compute_errors)
		{
			if (simVars.disc.timestepping_method.find("forcing")!=std::string::npos || simVars.disc.timestepping_method.find("_mms")!=std::string::npos)
			{
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
			}
			else //if (simVars.setup.benchmark_scenario_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppers.ln_cole_hopf->run_timestep(
						 ts_u, ts_v,
						 ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppers.l_direct->run_timestep(
						 ts_u, ts_v,
						 ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).reduce_rms_quad();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).reduce_rms_quad();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).reduce_maxAbs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).reduce_maxAbs();

		}
	}


	PlaneData compute_errors2(
         const PlaneData &i_planeData_u,
         const PlaneData &i_planeData_v
	)
	{
		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData u = i_planeData_u;
		PlaneData v = i_planeData_v;

		// Analytical solution at current time on original grid
		PlaneData ts_u = t0_prog_u;
		PlaneData ts_v = t0_prog_v;

		if (simVars.misc.compute_errors)
		{
			if (simVars.disc.timestepping_method.find("forcing")!=std::string::npos || simVars.disc.timestepping_method.find("_mms")!=std::string::npos)
			{
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
			}
			else //if (simVars.setup.benchmark_scenario_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppers.ln_cole_hopf->run_timestep(
						 ts_u, ts_v,
						 ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppers.l_direct->run_timestep(
						 ts_u, ts_v,
						 ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).reduce_rms_quad();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).reduce_rms_quad();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).reduce_maxAbs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).reduce_maxAbs();

			return ts_u;
		}
		return nullptr;
	}


	/*
	 * Check for final time step
	 */
	bool should_quit()
	{
		if (simVars.timecontrol.max_timesteps_nr != -1 && simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr)
			return true;

		if (!std::isinf(simVars.timecontrol.max_simulation_time))
			if (simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time+simVars.timecontrol.max_simulation_time*1e-10)	// care about roundoff errors with 1e-10
				return true;

		return false;
	}

#if SWEET_GUI

	/**
	 * Postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
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
		PlaneData ts_u = t0_prog_u;
		PlaneData ts_v = t0_prog_v;

		timeSteppers.ln_cole_hopf->run_timestep(
				ts_u, ts_v,
				ts_u, ts_v,
				simVars.timecontrol.current_simulation_time,
				0
		);

		switch(simVars.misc.vis_id)
		{
		case -1:
			vis = ts_u;
			break;

		case -2:
			vis = ts_u-prog_u;
			break;

		case -3:
			vis = ts_v-prog_v;
			break;
		}

		*o_dataArray = &vis;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];
		return;
	}


	/**
	 * Return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostic values if required
		update_diagnostics();

		const char* description = "";
		if (simVars.misc.vis_id >= 0)
		{
			int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
			description = vis_arrays[id].description;
		}
		else
		{
			switch (simVars.misc.vis_id)
			{
			case -1:
				description = "Cole-Hopf solution for u";
				break;

			case -2:
				description = "Diff in u to exact solution";
				break;

			case -3:
				description = "Diff in v to exact solution";
				break;
			}
		}

		static char title_string[1024];
		sprintf(title_string, "Time: %f, k: %i, dt: %.3e, Vis: %s, TEnergy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_energy,
				vis.reduce_max(),
				vis.reduce_min() );
		return title_string;
	}


	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;

		case 'c':
			// dump data arrays
			prog_u.file_physical_saveData_ascii("burgers_dump_u.csv");
			prog_v.file_physical_saveData_ascii("burgers_dump_v.csv");
			break;

		case 'C':
			// dump data arrays to VTK
			prog_u.file_physical_saveData_vtk("burgers_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("burgers_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_u.file_physical_loadData("burgers_dump_u.csv", simVars.setup.input_data_binary);
			prog_v.file_physical_loadData("burgers_dump_v.csv", simVars.setup.input_data_binary);
			break;
		}
	}
#endif


	bool instability_detected()
	{
		// Necessary to circumvent FFTW transformations on prog_u and prog_v, which would lead to errors
		PlaneData u = prog_u;
		PlaneData v = prog_v;

		return !(	u.reduce_boolean_all_finite() &&
					v.reduce_boolean_all_finite()
				);
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

//		timeSteppers.setup(
//				simVars.disc.timestepping_method,
//				simVars.disc.timestepping_order,
//				simVars.disc.timestepping_order2,
//				op,
//				simVars
//			);

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

		// Setup look-up table of Manufactured Solution
		if (simVars.disc.timestepping_method.find("_mms")!=std::string::npos)
		{
			timeSteppers.setup_look_up_table(
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.max_simulation_time,
				simVars.timecontrol.current_timestep_size
			);
			timeSteppersCoarse.setup_look_up_table(
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.max_simulation_time,
				simVars.timecontrol.current_timestep_size
			);
		}

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
		double dt = simVars.timecontrol.current_timestep_size;

		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// Reset dt as it might be changed due to timesize limitation
		simVars.timecontrol.current_timestep_size = dt;

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

		if (simVars.parareal.coarse_timestep_size <= 0)
			simVars.parareal.coarse_timestep_size = timeframe_end-timeframe_start;
		double tstep = simVars.parareal.coarse_timestep_size;

		// make multiple time steps in the coarse solver possible
		while (simVars.timecontrol.current_simulation_time < timeframe_end)
		{
			if (simVars.timecontrol.current_simulation_time + tstep > timeframe_end)
				tstep = timeframe_end - simVars.timecontrol.current_simulation_time;

			//run_timestep(simVars.parareal.coarse_timestepping_method, simVars.parareal.coarse_timestepping_order);
			timeSteppersCoarse.master->run_timestep(
				prog_u, prog_v,
				prog_u_prev, prog_v_prev,
				tstep,
				simVars.timecontrol.current_simulation_time
			);
			// Provide information to parameters
//			simVars.timecontrol.current_timestep_size = simVars.timecontrol.current_timestep_size;
			simVars.timecontrol.current_simulation_time += tstep;
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

		if (simVars.misc.compute_errors){
			compute_errors(prog_u, prog_v);
			std::cout << "maxabs error compared to analytical solution: " << benchmark.benchmark_analytical_error_maxabs_u << std::endl;
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


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_parareal(
			const PlaneData &i_planeData,
			const char* i_name, ///< name of output variable
			const int i_iteration_id,
			const int i_time_slice_id
		)
	{
		char buffer[1024];

		sprintf(buffer, (simVars.misc.output_file_name_prefix+std::string("output_%s_iter_%d_slice_%d.csv")).c_str(), i_name, i_iteration_id, i_time_slice_id);
		i_planeData.file_physical_saveData_ascii(buffer,'\n',12,1);

		/*
		char tmp[128];
		strcpy(tmp,i_name);
		strcat(tmp,"_spec");
		sprintf(buffer, (simVars.misc.output_file_name_prefix+std::string("output_%s_iter_%d_slice_%d.csv")).c_str(), tmp, simVars.timecontrol.current_simulation_time*simVars.misc.output_time_scale);
		i_planeData.file_spectral_saveData_ascii(buffer,'\n',12,1);
		*/

		return buffer;
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

		if (simVars.misc.output_file_name_prefix.size() > 0)
		{
			Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2>& data = (Parareal_Data_PlaneData<NUM_OF_UNKNOWNS*2>&)i_data;

			// Copy to minimize errors from FFT
			PlaneData tmp_u = *data.data_arrays[0];
			PlaneData tmp_v = *data.data_arrays[1];

			write_file_parareal(tmp_u,"prog_u",iteration_id,time_slice_id);

			char buffer[1024];
			sprintf(buffer,(simVars.misc.output_file_name_prefix+std::string("output_%s_iter_%d_slice_%d.csv")).c_str(), "prog_u_amp_phase", iteration_id, time_slice_id);

			std::ofstream file(buffer, std::ios_base::trunc);
			file << std::setprecision(12);

			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				file << x << ", " << tmp_u.spectral_return_amplitude(0,x) << ", " << tmp_u.spectral_return_phase(0,x) << std::endl;
			}
			file.close();
			file.clear();


			if (simVars.misc.compute_errors)
			{
				PlaneData ana = compute_errors2(tmp_u, tmp_v);

				write_file_parareal(ana,"analytical",iteration_id,time_slice_id);

				sprintf(buffer,(simVars.misc.output_file_name_prefix+std::string("output_%s_iter_%d_slice_%d.csv")).c_str(),"analytical_amp_phase", iteration_id, time_slice_id);

				file.open(buffer, std::ios_base::trunc);
				file << std::setprecision(12);
				for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
				{
					file << x << ", " << ana.spectral_return_amplitude(0,x) << ", " << ana.spectral_return_phase(0,x) << std::endl;
				}
				file.close();
			}
		}
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

		std::stringstream header;
		std::stringstream rows;

		rows << std::setprecision(12);

		// Prefix
		if (iteration_id == 0 && time_slice_id == 0)
			header << "DATA";
		rows << "DATA";

		// Time
		if (iteration_id == 0 && time_slice_id == 0)
			header << "\tITERATION\tTIME_SLICE";
		rows << "\t" << iteration_id << "\t" << time_slice_id;

		// Energy
		if (iteration_id == 0 && time_slice_id == 0)
			header << "\tTOTAL_ENERGY";
		rows << "\t" << simVars.diag.total_energy;

		if (simVars.misc.compute_errors)
		{
			// Copy to minimize errors from FFT
			PlaneData tmp_u = prog_u;
			if (iteration_id == 0 && time_slice_id == 0)
				header << "\tMAX_ABS_U\tMAX_RMS_U\tMAX_U";
			rows << "\t" << benchmark.benchmark_analytical_error_maxabs_u << "\t" << benchmark.benchmark_analytical_error_rms_u << "\t" << tmp_u.reduce_max();
		}

		if (iteration_id == 0 && time_slice_id == 0)
			std::cout << header.str() << std::endl;

		std::cout << rows.str() << std::endl;

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
			nullptr
	};

	// default values for program specific parameter

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << std::endl;

/*
 * #if SWEET_PARAREAL
 * 		simVars.parareal.printOptions();
 * #endif
 */
		return -1;
	}

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	// Print header
	std::cout << "---------------------------------" << std::endl;
	std::cout << "Solving viscous Burgers' equation" << std::endl;
	std::cout << "---------------------------------" << std::endl;
	simVars.outputConfig();
	std::cout << "LOCAL PARAMETERS" << std::endl;
	std::cout << "Computing error: " << simVars.misc.compute_errors << std::endl;
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
			// already called in constructor
//			simulationBurgers->reset();

			// Time counter
			Stopwatch time;

			// Start counting time
			time.reset();

			// Main time loop
			while(true)
			{
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
				std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs(simVars.diag.total_energy-simulationBurgers->diagnostics_energy_start) << std::endl;

				if (simVars.misc.compute_errors)
				{
					std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationBurgers->benchmark.benchmark_analytical_error_maxabs_u << std::endl;
					std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationBurgers->benchmark.benchmark_analytical_error_maxabs_v << std::endl;
				}
			}

			delete simulationBurgers;
		}
	}

	return 0;
}
