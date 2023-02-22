/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/burgers_timeintegrators
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/burgers_benchmarks
 */

#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif

#include <sweet/Stopwatch.hpp>
#include <sweet/SWEETError.hpp>


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData_Physical.hpp>
#include <sweet/plane/PlaneData_Spectral.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/Convert_ScalarDataArray_to_PlaneDataPhysical.hpp>

#include "burgers_timeintegrators/Burgers_Plane_TimeSteppers.hpp"

#include "burgers_benchmarks/BurgersValidationBenchmarks.hpp"




#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#if SWEET_PARAREAL
#	include <parareal/Parareal.hpp>
#endif


const int NUM_OF_UNKNOWNS=2;

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

// General parameters
SimulationVariables simVars;



class SimulationInstance
////#if SWEET_PARAREAL
////		:
////		public Parareal_SimulationInstance
////#endif
{

public:
	// Prognostic variables
	// u: velocity in x-direction
	// v: velocity in y-direction
	PlaneData_Spectral prog_u, prog_v;

	// Prognostic variables at time step t-dt
	PlaneData_Spectral prog_u_prev, prog_v_prev;

#if SWEET_GUI
	// visualization variable
	PlaneData_Physical vis;
#endif

	// Initial values for comparison with analytical solution
	PlaneData_Spectral t0_prog_u, t0_prog_v;

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
		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		// Calls initialization of the run (e.g. sets u, v)
		reset();

	}


	virtual ~SimulationInstance()
	{
	}



	void reset()
	{
#if 0
		if (simVars.benchmark.benchmark_id <0)
		{
			std::cout << std::endl;
			std::cout << "Benchmark scenario not selected (option -s [id])" << std::endl;
			BurgersValidationBenchmarks::printScenarioInformation();
			SWEETError("Benchmark scenario not selected");
		}
#endif
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
//#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_u.spectral_set_zero();
		prog_v.spectral_set_zero();
		///prog_u_prev.spectral_set_zero();
		///prog_v_prev.spectral_set_zero();
////#endif
////		prog_u.physical_set_all(0);
////		prog_v.physical_set_all(0);
////		prog_u_prev.physical_set_all(0);
////		prog_v_prev.physical_set_all(0);

		//Check if input parameters are adequate for this simulation
		if (simVars.disc.space_grid_use_c_staggering && simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.space_grid_use_c_staggering ||  !simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
#endif

		PlaneData_Physical u_phys(planeDataConfig);
		PlaneData_Physical v_phys(planeDataConfig);

		// Set initial conditions given from BurgersValidationBenchmarks
		if (simVars.disc.space_grid_use_c_staggering)
		{
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);
			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
				io_data = 0.0;
#if 0
					double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
					double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
				}
			);
		}
		else
		{
			u_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
					double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
				}
			);

			v_phys.physical_update_lambda_array_indices(
				[&](int i, int j, double &io_data)
				{
				io_data = 0.0;
#if 0
					double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
					double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
					io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
				}
			);
		}

		prog_u.loadPlaneDataPhysical(u_phys);
		prog_v.loadPlaneDataPhysical(v_phys);

		// Initialize t-dt time step with initial condition
		//prog_u_prev = prog_u;
		//prog_v_prev = prog_v;
		// Save initial conditions for analytic solution
		t0_prog_u = prog_u;
		t0_prog_v = prog_v;

		// Load data, if requested
		if (simVars.iodata.initial_condition_data_filenames.size() > 0)
			prog_u.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[0].c_str(), simVars.iodata.initial_condition_input_data_binary);

		if (simVars.iodata.initial_condition_data_filenames.size() > 1)
			prog_v.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[1].c_str(), simVars.iodata.initial_condition_input_data_binary);


		timeSteppers.setup(
			simVars.disc.timestepping_method,
			simVars.disc.timestepping_order,
			simVars.disc.timestepping_order2,
			op,
			simVars
		);

		update_diagnostics();
		diagnostics_energy_start = simVars.diag.total_energy;

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
            SWEETError("Computing errors for this timestepping-method is not possible");
      }
	}


	// Calculate the model diagnostics
	void update_diagnostics()
	{
		// Assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;

		double normalization = (simVars.sim.plane_domain_size[0]*simVars.sim.plane_domain_size[1]) /
								((double)simVars.disc.space_res_physical[0]*(double)simVars.disc.space_res_physical[1]);

		// Reduce amount of possible FFTs to minimize numerical error
		PlaneData_Physical tmp_u = prog_u.toPhys();
		PlaneData_Physical tmp_v = prog_v.toPhys();

		// Energy
		simVars.diag.total_energy =
			0.5*((
					tmp_u*tmp_u +
					tmp_v*tmp_v
				).physical_reduce_sum_quad()) * normalization;
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
				///prog_u_prev, prog_v_prev,
				simVars.timecontrol.current_timestep_size+simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time
			);
		}
		else
		{
			timeSteppers.master->run_timestep(
				prog_u, prog_v,
				///prog_u_prev, prog_v_prev,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);
		}

		// Advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			SWEETError("Max simulation time exceeded!");

		timestep_output();
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const PlaneData_Physical &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		char buffer[1024];

		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		//i_planeData.file_physical_saveData_ascii(buffer, '\n', 12, 1);
		i_planeData.file_physical_saveData_ascii(buffer);

		return buffer;
	}



	std::string output_filenames;

public:
	bool timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
		// output each time step
		if (simVars.iodata.output_each_sim_seconds < 0)
			return false;

		if (simVars.iodata.output_next_sim_seconds-simVars.iodata.output_next_sim_seconds*(1e-12) > simVars.timecontrol.current_simulation_time)
			return false;

		PlaneData_Physical u_phys = prog_u.toPhys();
		PlaneData_Physical v_phys = prog_v.toPhys();

		// Dump data in csv, if requested
		if (simVars.iodata.output_file_name.size() > 0)
		{
			output_filenames = "";
			output_filenames += write_file(u_phys, "prog_u");
			write_file(v_phys, "prog_v");
			//output_filenames += write_file(v_phys, "prog_v");
			//write_file(tmp_v, "prog_v");

			char buffer[1024];
			const char* filename_template = simVars.iodata.output_file_name.c_str();
			sprintf(buffer,filename_template,"prog_u_amp_phase",simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);

			std::ofstream file(buffer, std::ios_base::trunc);
			file << std::setprecision(12);

			for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
			{
				file << x << ", " << prog_u.spectral_return_amplitude(0,x) << ", " << prog_u.spectral_return_phase(0,x) << std::endl;
			}
			file.close();
		}

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();
			if (simVars.misc.compute_errors)
			{
				PlaneData_Physical tmp(planeDataConfig);
				tmp = compute_errors2(prog_u, prog_v).toPhys();

				write_file(tmp, "analytical");

				PlaneData_Spectral tmp_spec(planeDataConfig);
				tmp_spec.loadPlaneDataPhysical(tmp);

				char buffer[1024];
				const char* filename_template = simVars.iodata.output_file_name.c_str();
				sprintf(buffer,filename_template,"analytical_amp_phase",simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);

				std::ofstream file(buffer, std::ios_base::trunc);
				file << std::setprecision(12);
				for (std::size_t x = 0; x < planeDataConfig->spectral_data_size[0]; x++)
				{
					file << x << ", " << tmp_spec.spectral_return_amplitude(0,x) << ", " << tmp_spec.spectral_return_phase(0,x) << std::endl;
				}
				file.close();
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
				rows << "\t" << benchmark.benchmark_analytical_error_maxabs_u << "\t" << benchmark.benchmark_analytical_error_rms_u << "\t" << prog_u.spectral_reduce_max_abs();
			}

			if (simVars.timecontrol.current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;
		}

		if (simVars.iodata.output_each_sim_seconds > 0)
		{
			if (simVars.iodata.output_next_sim_seconds == simVars.timecontrol.max_simulation_time)
			{
				simVars.iodata.output_next_sim_seconds = std::numeric_limits<double>::infinity();
			}
			else
			{
				while (simVars.iodata.output_next_sim_seconds <= simVars.timecontrol.current_simulation_time)
					simVars.iodata.output_next_sim_seconds += simVars.iodata.output_each_sim_seconds;

				if (simVars.iodata.output_next_sim_seconds > simVars.timecontrol.max_simulation_time)
					simVars.iodata.output_next_sim_seconds = simVars.timecontrol.max_simulation_time;
			}
		}

		return true;
	}


public:
	void compute_errors(
         const PlaneData_Spectral &i_planeData_u,
         const PlaneData_Spectral &i_planeData_v
	)
	{
		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData_Physical u = i_planeData_u.toPhys();
		PlaneData_Physical v = i_planeData_v.toPhys();

		// Analytical solution at current time on original grid
		PlaneData_Spectral ts_u = t0_prog_u;
		PlaneData_Spectral ts_v = t0_prog_v;

		PlaneData_Physical ts_u_phys = ts_u.toPhys();
		PlaneData_Physical ts_v_phys = ts_v.toPhys();

		if (simVars.misc.compute_errors)
		{
			//if (simVars.setup.benchmark_id > 51 && simVars.setup.benchmark_id < 65)
			if (simVars.disc.timestepping_method.find("forcing")!=std::string::npos)
			{
				if (simVars.disc.space_grid_use_c_staggering)
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				else
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				ts_u.loadPlaneDataPhysical(ts_u_phys);
				ts_v.loadPlaneDataPhysical(ts_v_phys);
			}
			else //if (simVars.setup.benchmark_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppers.ln_cole_hopf->run_timestep(
						 ts_u, ts_v,
						 ////ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppers.l_direct->run_timestep(
						 ts_u, ts_v,
						 ////ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).toPhys().physical_reduce_rms();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).toPhys().physical_reduce_rms();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).toPhys().physical_reduce_max_abs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).toPhys().physical_reduce_max_abs();

		}
	}


	PlaneData_Spectral compute_errors2(
         const PlaneData_Spectral &i_planeData_u,
         const PlaneData_Spectral &i_planeData_v
	)
	{
		// Necessary to circumvent FFTW transformations on i_planeData_u and i_planeData_v, which would lead to errors
		PlaneData_Physical u = i_planeData_u.toPhys();
		PlaneData_Physical v = i_planeData_v.toPhys();

		// Analytical solution at current time on original grid
		PlaneData_Spectral ts_u = t0_prog_u;
		PlaneData_Spectral ts_v = t0_prog_v;

		PlaneData_Physical ts_u_phys = ts_u.toPhys();
		PlaneData_Physical ts_v_phys = ts_v.toPhys();

		if (simVars.misc.compute_errors)
		{
			//if (simVars.setup.benchmark_id > 51 && simVars.setup.benchmark_id < 65)
			if (simVars.disc.timestepping_method.find("forcing")!=std::string::npos)
			{
				if (simVars.disc.space_grid_use_c_staggering)
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];
							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				else
				{
					ts_u_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_u(simVars, x, y);
						}
					);

					ts_v_phys.physical_update_lambda_array_indices(
						[&](int i, int j, double &io_data)
						{
							io_data = 0.0;
#if 0
							double x = (((double)i+0.5)/(double)simVars.disc.space_res_physical[0])*simVars.sim.plane_domain_size[0];
							double y = (((double)j+0.5)/(double)simVars.disc.space_res_physical[1])*simVars.sim.plane_domain_size[1];

							io_data = BurgersValidationBenchmarks::return_v(simVars, x, y);
#endif
						}
					);
				}
				ts_u.loadPlaneDataPhysical(ts_u_phys);
				ts_v.loadPlaneDataPhysical(ts_v_phys);
			}
			else //if (simVars.setup.benchmark_id == 70)
			{
				if (analytic_solution == 1)
				{
				   timeSteppers.ln_cole_hopf->run_timestep(
						 ts_u, ts_v,
						 /////ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
				else if (analytic_solution == 2)
				{
				   timeSteppers.l_direct->run_timestep(
						 ts_u, ts_v,
						 /////ts_u, ts_v,
						 simVars.timecontrol.current_simulation_time,
						 0
				   );
				}
			}
			benchmark.benchmark_analytical_error_rms_u = (ts_u-u).toPhys().physical_reduce_rms();
			benchmark.benchmark_analytical_error_rms_v = (ts_v-v).toPhys().physical_reduce_rms();

			benchmark.benchmark_analytical_error_maxabs_u = (ts_u-u).toPhys().physical_reduce_max_abs();
			benchmark.benchmark_analytical_error_maxabs_v = (ts_v-v).toPhys().physical_reduce_max_abs();

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
		const PlaneData_Physical* data;
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
			const PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		PlaneData_Physical ts_u = t0_prog_u.toPhys();
		PlaneData_Physical ts_v = t0_prog_v.toPhys();

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
		*o_aspect_ratio = simVars.sim.plane_domain_size[1] / simVars.sim.plane_domain_size[0];
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
				vis.physical_reduce_max(),
				vis.physical_reduce_min() );
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

#if 0
		case 'l':
			// load data arrays
			prog_u.file_physical_loadData("burgers_dump_u.csv", simVars.benchmark.initial_condition_input_data_binary);
			prog_v.file_physical_loadData("burgers_dump_v.csv", simVars.benchmark.initial_condition_input_data_binary);
#endif

			break;
		}
	}
#endif


	bool instability_detected()
	{
		// Necessary to circumvent FFTW transformations on prog_u and prog_v, which would lead to errors
		PlaneData_Spectral u = prog_u;
		PlaneData_Spectral v = prog_v;

		return !(	u.toPhys().physical_reduce_boolean_all_finite() &&
					v.toPhys().physical_reduce_boolean_all_finite()
				);
	}


};




int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif
	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << std::endl;

		return -1;
	}

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

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

			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			// Set planeDataConfig and planeOperators for each level
			std::vector<PlaneDataConfig*> planeDataConfigs;
			std::vector<PlaneOperators*> ops;

			// fine
			planeDataConfigs.push_back(planeDataConfig);
			ops.push_back(&op);

			// coarse
			if (simVars.parareal.spatial_coarsening)
			{
				///for (int j = 0; j < 2; j++)
				///	assert(simVars.disc.space_res_physical[j] == -1);
				int N_physical[2] = {-1, -1};
				int N_spectral[2];
				double frac;
				if ( simVars.parareal.coarse_timestep_size > 0)
					frac = simVars.timecontrol.current_timestep_size / simVars.parareal.coarse_timestep_size;
				else
					frac = simVars.timecontrol.current_timestep_size / (simVars.timecontrol.max_simulation_time / simVars.parareal.coarse_slices );
				for (int j = 0; j < 2; j++)
					N_spectral[j] = std::max(4, int(simVars.disc.space_res_spectral[j] * frac));
				planeDataConfigs.push_back(new PlaneDataConfig);
				planeDataConfigs.back()->setupAuto(N_physical, N_spectral, simVars.misc.reuse_spectral_transformation_plans);

				ops.push_back(new PlaneOperators(planeDataConfigs.back(), simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs));
			}
			else
			{
				planeDataConfigs.push_back(planeDataConfig);
				ops.push_back(&op);
			}

			Burgers_Plane_TimeSteppers* timeSteppersFine = new Burgers_Plane_TimeSteppers;
			Burgers_Plane_TimeSteppers* timeSteppersCoarse = new Burgers_Plane_TimeSteppers;

			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller<Burgers_Plane_TimeSteppers, 2> parareal_Controller(&simVars,
												planeDataConfigs,
												ops,
												timeSteppersFine,
												timeSteppersCoarse);

			// setup controller. This initializes several simulation instances
			parareal_Controller.setup();

			// execute the simulation
			parareal_Controller.run();

			delete timeSteppersFine;
			delete timeSteppersCoarse;

			if (simVars.parareal.spatial_coarsening)
			{
				delete planeDataConfigs[1];
				delete ops[1];
			}
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

			// Start counting time
			time.start();

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

			// Stop counting time
			time.stop();

			double seconds = time();

			if (simVars.iodata.output_file_name.size() > 0)
				std::cout << "[MULE] reference_filenames: " << simulationBurgers->output_filenames << std::endl;

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
