/*
 * SWE with nonlinear part using REXI test
 *
 */

#if SWEET_GUI
	#include <sweet/VisSweet.hpp>
#endif

#ifndef SWEET_MPI
	// just for programming purpose
	#define  SWEET_MPI 1
#endif


#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneData.hpp>

#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDataSampler.hpp>
#include <sweet/plane/PlaneDataGridMapping.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include <sweet/plane/Convert_PlaneDataComplex_to_PlaneData.hpp>
#include <sweet/plane/Convert_PlaneData_to_PlaneDataComplex.hpp>
#include <sweet/Stopwatch.hpp>
#include <sweet/FatalError.hpp>
#include <benchmarks_plane/SWEPlaneBenchmarksCombined.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include <sweet/SimulationBenchmarkTiming.hpp>


#include "swe_plane/SWE_Plane_TimeSteppers.hpp"
#include "swe_plane/SWE_Plane_Normal_Modes.hpp"


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;

/// general parameters
SimulationVariables simVars;

double param_initial_freq_x_mul;
double param_initial_freq_y_mul;


class SimulationInstance
{
public:
	// Prognostic variables
	// h: surface height (perturbation)
	// u: velocity in x-direction
	// v: velocity in y-direction
	PlaneData prog_h_pert, prog_u, prog_v;


	// Initial values for comparison with analytical solution
	PlaneData t0_prog_h_pert, t0_prog_u, t0_prog_v;

	// Forcings
	PlaneData force_h_pert, force_u, force_v;

	PlaneDataGridMapping gridMapping;


	// implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppers;

	SWEPlaneBenchmarksCombined swePlaneBenchmarks;

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	class BenchmarkErrors
	{
	public:
		//Max difference to initial conditions
		double t0_error_max_abs_h_pert;
		double t0_error_max_abs_u;
		double t0_error_max_abs_v;

		// Error measures L2 norm
		double analytical_error_rms_h;
		double analytical_error_rms_u;
		double analytical_error_rms_v;

		// Error measures max norm
		double analytical_error_maxabs_h;
		double analytical_error_maxabs_u;
		double analytical_error_maxabs_v;
	};

	BenchmarkErrors benchmark;

	// Finite difference operators
	PlaneOperators op;


	/// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;
	double diagnostics_mass_start = 0;
	double diagnostics_potential_entrophy_start = 0;


	bool compute_error_difference_to_initial_condition = false;
	bool compute_error_to_analytical_solution = false;



public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h_pert(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

		t0_prog_h_pert(planeDataConfig),
		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

		force_h_pert(planeDataConfig),
		force_u(planeDataConfig),
		force_v(planeDataConfig),

		// Initialises operators
		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
	{
		// Calls initialisation of the run (e.g. sets u, v, h)
		reset();
	}


	virtual ~SimulationInstance()
	{
	}



	void reset()
	{
		SimulationBenchmarkTimings::getInstance().main_setup.start();

		simVars.reset();

		if (simVars.benchmark.benchmark_name == "")
		{
			std::cout << "Benchmark scenario not selected (option --benchmark-name [string])" << std::endl;
			swePlaneBenchmarks.printBenchmarkInformation();
			FatalError("Benchmark name not given");
		}


		// Initialise diagnostics
		last_timestep_nr_update_diagnostics = -1;


		benchmark.t0_error_max_abs_h_pert = -1;
		benchmark.t0_error_max_abs_u = -1;
		benchmark.t0_error_max_abs_v = -1;

		benchmark.analytical_error_rms_h = -1;
		benchmark.analytical_error_rms_u = -1;
		benchmark.analytical_error_rms_v = -1;

		benchmark.analytical_error_maxabs_h = -1;
		benchmark.analytical_error_maxabs_u = -1;
		benchmark.analytical_error_maxabs_v = -1;

		simVars.timecontrol.current_timestep_nr = 0;
		simVars.timecontrol.current_simulation_time = 0;

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_h_pert.spectral_set_all(0, 0);
		prog_u.spectral_set_all(0, 0);
		prog_v.spectral_set_all(0, 0);
#endif

		// Setup prog vars
		prog_h_pert.physical_set_all(simVars.sim.h0);
		prog_u.physical_set_all(0);
		prog_v.physical_set_all(0);

		// Check if input parameters are adequate for this simulation
		if (simVars.disc.space_grid_use_c_staggering && simVars.disc.space_use_spectral_basis_diffs)
			FatalError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.space_grid_use_c_staggering ||  !simVars.disc.space_use_spectral_basis_diffs)
			FatalError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
#endif


		if (simVars.disc.space_grid_use_c_staggering)
			gridMapping.setup(simVars, planeDataConfig);

		swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, simVars, op);

		prog_h_pert = t0_prog_h_pert;
		prog_u = t0_prog_u;
		prog_v = t0_prog_v;


		// Load data, if requested
		if (simVars.iodata.initial_condition_data_filenames.size() > 0)
			prog_h_pert.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[0].c_str(), simVars.iodata.initial_condition_input_data_binary);

		if (simVars.iodata.initial_condition_data_filenames.size() > 1)
			prog_u.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[1].c_str(), simVars.iodata.initial_condition_input_data_binary);

		if (simVars.iodata.initial_condition_data_filenames.size() > 2)
			prog_v.file_physical_loadData(simVars.iodata.initial_condition_data_filenames[2].c_str(), simVars.iodata.initial_condition_input_data_binary);

		timeSteppers.setup(
				simVars.disc.timestepping_method,
				simVars.disc.timestepping_order,
				simVars.disc.timestepping_order2,
				op,
				simVars
			);

		if (simVars.misc.compute_errors)
		{
			//Compute difference to initial condition (makes more sense in steady state cases, but useful in others too)
			compute_error_difference_to_initial_condition = true;
					//simVars.setup.benchmark_id == 2 ||
					//simVars.setup.benchmark_id == 3 ||
					//simVars.setup.benchmark_id == 14;

			//Compute difference to analytical solution (makes more sense in linear cases, but might be useful in others too)
			compute_error_to_analytical_solution = timeSteppers.linear_only;
		}
		else
		{
			compute_error_difference_to_initial_condition = false;
			compute_error_to_analytical_solution = false;
		}

		update_diagnostics();

		diagnostics_energy_start = simVars.diag.total_energy;
		diagnostics_mass_start = simVars.diag.total_mass;
		diagnostics_potential_entrophy_start = simVars.diag.total_potential_enstrophy;


		timestep_do_output();

		SimulationBenchmarkTimings::getInstance().main_setup.stop();

	}


	//Calculate the model diagnostics
	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == simVars.timecontrol.current_timestep_nr)
			return;

		last_timestep_nr_update_diagnostics = simVars.timecontrol.current_timestep_nr;


		if (simVars.disc.space_grid_use_c_staggering)
		{
			PlaneDiagnostics::update_staggered_huv_to_mass_energy_enstrophy(
					op,
					prog_h_pert,
					prog_u,
					prog_v,
					simVars
			);
		}
		else
		{
			PlaneDiagnostics::update_nonstaggered_huv_to_mass_energy_enstrophy(
					op,
					prog_h_pert,
					prog_u,
					prog_v,
					simVars
			);
		}
	}



	void normal_mode_analysis()
	{
		SWE_Plane_Normal_Modes::normal_mode_analysis(
								prog_h_pert,
								prog_u,
								prog_v,
								3,
								simVars,
								this,
								&SimulationInstance::run_timestep
						);

		std::cout << "\n Done normal mode analysis in separate class" << std::endl;
	}


	/**
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		if (simVars.timecontrol.current_simulation_time + simVars.timecontrol.current_timestep_size > simVars.timecontrol.max_simulation_time)
			simVars.timecontrol.current_timestep_size = simVars.timecontrol.max_simulation_time - simVars.timecontrol.current_simulation_time;

		timeSteppers.master->run_timestep(
				prog_h_pert, prog_u, prog_v,
				simVars.timecontrol.current_timestep_size,
				simVars.timecontrol.current_simulation_time
			);

		// Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral

		if (simVars.sim.viscosity != 0 && simVars.misc.use_nonlinear_only_visc == 0)
		{
#if !SWEET_USE_PLANE_SPECTRAL_SPACE //TODO: this needs checking

			double dt = simVars.timecontrol.current_timestep_size;

			prog_u = prog_u + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_u, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_v = prog_v + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_v, simVars.sim.viscosity_order)*simVars.sim.viscosity;
			prog_h_pert = prog_h_pert + pow(-1,simVars.sim.viscosity_order/2)* dt*op.diffN_x(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity
					+ pow(-1,simVars.sim.viscosity_order/2)*dt*op.diffN_y(prog_h_pert, simVars.sim.viscosity_order)*simVars.sim.viscosity;
#else
			prog_u = op.implicit_diffusion(prog_u, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_v = op.implicit_diffusion(prog_v, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
			prog_h_pert = op.implicit_diffusion(prog_h_pert, simVars.timecontrol.current_timestep_size*simVars.sim.viscosity, simVars.sim.viscosity_order);
#endif
		}

		// advance time step and provide information to parameters
		simVars.timecontrol.current_simulation_time += simVars.timecontrol.current_timestep_size;
		simVars.timecontrol.current_timestep_nr++;

		if (simVars.timecontrol.current_simulation_time > simVars.timecontrol.max_simulation_time)
			FatalError("Max simulation time exceeded!");

		timestep_do_output();
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

		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write spectrum info to data and return string of file name
	 */

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	std::string write_file_spec(
			const PlaneData &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		char buffer[1024];

		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		i_planeData.file_spectral_abs_saveData_ascii(buffer);
		//i_planeData.file_spectral_saveData_ascii(buffer);
		return buffer;
	}
#endif


	std::string output_filenames;


public:
	bool timestep_do_output(
			std::ostream &o_ostream = std::cout
	)
	{
		if (simVars.misc.normal_mode_analysis_generation > 0)
			return false;

		// output each time step
		if (simVars.iodata.output_each_sim_seconds < 0)
			return false;

		if (simVars.iodata.output_next_sim_seconds-simVars.iodata.output_next_sim_seconds*(1e-12) > simVars.timecontrol.current_simulation_time)
			return false;

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		{
			// For output, variables need to be on unstaggered A-grid
			PlaneData t_h(planeDataConfig);
			PlaneData t_u(planeDataConfig);
			PlaneData t_v(planeDataConfig);

			if (simVars.disc.space_grid_use_c_staggering) // Remap in case of C-grid
			{
				t_h = prog_h_pert;
				gridMapping.mapCtoA_u(prog_u, t_u);
				gridMapping.mapCtoA_v(prog_v, t_v);
			}
			else
			{
				t_h = prog_h_pert;
				t_u = prog_u;
				t_v = prog_v;
			}

			//std::cout << simVars.inputoutput.output_next_sim_seconds << "\t" << simVars.timecontrol.current_simulation_time << std::endl;

			// Dump  data in csv, if output filename is not empty
			if (simVars.iodata.output_file_name.size() > 0)
			{
				output_filenames = "";

				output_filenames = write_file(t_h, "prog_h_pert");
				output_filenames += ";" + write_file(t_u, "prog_u");
				output_filenames += ";" + write_file(t_v, "prog_v");

				output_filenames += ";" + write_file(op.ke(t_u,t_v),"diag_ke");

#if SWEET_USE_PLANE_SPECTRAL_SPACE
				output_filenames += ";" + write_file_spec(op.ke(t_u,t_v),"diag_ke_spec");
#endif

				output_filenames += ";" + write_file(op.vort(t_u, t_v), "diag_vort");
				output_filenames += ";" + write_file(op.div(t_u, t_v), "diag_div");

			}
		}


		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();
			compute_errors();


			if (simVars.timecontrol.current_timestep_nr == 0)
			{
				simVars.diag.outputConfig();
			}

			std::stringstream header;
			std::stringstream rows;

			rows << std::setprecision(16);

			// Prefix
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "DATA";
			rows << "DATA";

			// Time
			if (simVars.timecontrol.current_timestep_nr == 0)
				header << "\tT";
			rows << "\t" << simVars.timecontrol.current_simulation_time;

#if 1
			// Mass, Energy, Enstrophy
			header << "\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";
			rows << "\t" << simVars.diag.total_mass;
			rows << "\t" << simVars.diag.total_energy;
			rows << "\t" << simVars.diag.total_potential_enstrophy;

			// Mass, Energy, Enstrophy
			header << "\tTOTAL_MASS_REL_ERROR\tTOTAL_ENERGY_REL_ERROR\tPOT_ENSTROPHY_REL_ERROR";
			rows << "\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start);
			rows << "\t" << std::abs((simVars.diag.total_energy-diagnostics_energy_start)/diagnostics_energy_start) ;
			rows << "\t" << std::abs((simVars.diag.total_potential_enstrophy-diagnostics_potential_entrophy_start)/diagnostics_potential_entrophy_start);
#endif

#if 1
			if (compute_error_difference_to_initial_condition)
			{
				// Difference to initial condition
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tDIFF_MAXABS_H0\tDIFF_MAXABS_U0\tDIFF_MAXABS_V0";

				rows << "\t" << benchmark.t0_error_max_abs_h_pert << "\t" << benchmark.t0_error_max_abs_u << "\t" << benchmark.t0_error_max_abs_v;
			}
#endif

#if 1
			if (compute_error_to_analytical_solution)
			{
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

				rows << "\t" << benchmark.analytical_error_maxabs_h << "\t" << benchmark.analytical_error_maxabs_u << "\t" << benchmark.analytical_error_maxabs_v;
			}
#endif

			if (simVars.timecontrol.current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;

#if 1
			if (diagnostics_mass_start > 0.00001 && std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) > 10000000.0)
			{
				std::cerr << "\n DIAGNOSTICS MASS DIFF TOO LARGE:\t" << std::abs((simVars.diag.total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			}
#endif

		}

		if (simVars.iodata.output_next_sim_seconds == simVars.timecontrol.max_simulation_time)
		{
			simVars.iodata.output_next_sim_seconds = std::numeric_limits<double>::infinity();
		}
		else
		{
			while (simVars.iodata.output_next_sim_seconds-simVars.iodata.output_next_sim_seconds*(1e-12) <= simVars.timecontrol.current_simulation_time)
				simVars.iodata.output_next_sim_seconds += simVars.iodata.output_each_sim_seconds;

			if (simVars.iodata.output_next_sim_seconds > simVars.timecontrol.max_simulation_time)
				simVars.iodata.output_next_sim_seconds = simVars.timecontrol.max_simulation_time;
		}

		return true;
	}


public:
	void compute_errors()
	{
		if (compute_error_difference_to_initial_condition || compute_error_to_analytical_solution)
		{
			/**
			 * Compute difference to initial condition
			 */
			if (compute_error_difference_to_initial_condition)
			{
				benchmark.t0_error_max_abs_h_pert = (prog_h_pert - t0_prog_h_pert).reduce_maxAbs();
				benchmark.t0_error_max_abs_u = (prog_u - t0_prog_u).reduce_maxAbs();
				benchmark.t0_error_max_abs_v = (prog_v - t0_prog_v).reduce_maxAbs();
			}

			// Calculate linear exact solution, if compute error requests
			if (compute_error_to_analytical_solution)
			{
				// Analytical solution at specific time on A-grid
				PlaneData ts_h_pert = t0_prog_h_pert;
				PlaneData ts_u = t0_prog_u;
				PlaneData ts_v = t0_prog_v;

				// Run exact solution for linear case
				timeSteppers.l_direct->run_timestep(
						ts_h_pert, ts_u, ts_v,
						simVars.timecontrol.current_simulation_time,	// time step size
						0				// initial condition given at time 0
				);

				benchmark.analytical_error_rms_h = (ts_h_pert-prog_h_pert).reduce_rms_quad();
				benchmark.analytical_error_rms_u = (ts_u-prog_u).reduce_rms_quad();
				benchmark.analytical_error_rms_v = (ts_v-prog_v).reduce_rms_quad();

				benchmark.analytical_error_maxabs_h = (ts_h_pert-prog_h_pert).reduce_maxAbs();
				benchmark.analytical_error_maxabs_u = (ts_u-prog_u).reduce_maxAbs();
				benchmark.analytical_error_maxabs_v = (ts_v-prog_v).reduce_maxAbs();
			}
		}
	}



public:
	bool should_quit()
	{
		if (
				simVars.timecontrol.max_timesteps_nr != -1 &&
				simVars.timecontrol.max_timesteps_nr <= simVars.timecontrol.current_timestep_nr
		)
			return true;

		if (!std::isinf(simVars.timecontrol.max_simulation_time))
			if (simVars.timecontrol.max_simulation_time <= simVars.timecontrol.current_simulation_time+simVars.timecontrol.max_simulation_time*1e-10)	// care about roundoff errors with 1e-10
				return true;

		return false;
	}


	bool instability_detected()
	{
		return !(	prog_h_pert.reduce_boolean_all_finite() &&
					prog_u.reduce_boolean_all_finite() &&
					prog_v.reduce_boolean_all_finite()
				);
	}


};



int main(int i_argc, char *i_argv[])
{

	MemBlockAlloc::setup();

	//input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"initial-freq-x-mul",		/// frequency multipliers for special scenario setup
			"initial-freq-y-mul",
			nullptr
	};

	// default values for specific input (for general input see SimulationVariables.hpp)
	//simVars.bogus.var[0];  //frequency in x for waves test case
	//simVars.bogus.var[1];  //frequency in y for waves test case

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << std::endl;
		std::cout << "Special parameters:" << std::endl;
		std::cout << "	--initial-freq-x-mul [float]" << std::endl;
		std::cout << "	--initial-freq-y-mul [float]" << std::endl;
		std::cout << "" << std::endl;

		return -1;
	}
	// Frequency for certain initial conditions
	param_initial_freq_x_mul = 0;
	param_initial_freq_y_mul = 0;

	if (simVars.bogus.var[0] != "")
		param_initial_freq_x_mul = atof(simVars.bogus.var[0].c_str());

	if (simVars.bogus.var[1] != "")
		param_initial_freq_y_mul = atof(simVars.bogus.var[1].c_str());

	if (simVars.misc.verbosity > 5)
		std::cout << " + Setting up FFT plans..." << std::flush;

	planeDataConfigInstance.setupAuto(simVars.disc.space_res_physical, simVars.disc.space_res_spectral, simVars.misc.reuse_spectral_transformation_plans);

	if (simVars.misc.verbosity > 5)
		std::cout << " done" << std::endl;

	planeDataConfig->printInformation();

	// Print header
	std::cout << std::endl;
	simVars.outputConfig();
	std::cout << "Computing error: " << simVars.misc.compute_errors << std::endl;
	std::cout << std::endl;

	std::ostringstream buf;
	buf << std::setprecision(14);
	{


		SimulationBenchmarkTimings::getInstance().main.start();

		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			// also initializes diagnostics
			// already called in constructor
			//simulationSWE->reset();



			if (simVars.misc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normal_mode_analysis();
			}
			else
			{
				SimulationBenchmarkTimings::getInstance().main_timestepping.start();

				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Main call for timestep run
					simulationSWE->run_timestep();

					// Instability
					if (simVars.misc.instability_checks)
					{
						if (simulationSWE->instability_detected())
							FatalError("INSTABILITY DETECTED");
					}
				}

				SimulationBenchmarkTimings::getInstance().main_timestepping.stop();
			}


			if (simVars.iodata.output_file_name.size() > 0)
				std::cout << "[MULE] reference_filenames: " << simulationSWE->output_filenames << std::endl;

			// End of run output results
			std::cout << "***************************************************" << std::endl;
			std::cout << "Number of time steps: " << simVars.timecontrol.current_timestep_nr << std::endl;
			std::cout << "Time per time step: " << SimulationBenchmarkTimings::getInstance().main_timestepping()/(double)simVars.timecontrol.current_timestep_nr << " sec/ts" << std::endl;
			std::cout << "Last time step size: " << simVars.timecontrol.current_timestep_size << std::endl;

			simulationSWE->compute_errors();


			{
				if (simVars.misc.verbosity > 0)
				{
					std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-simulationSWE->diagnostics_energy_start)/simulationSWE->diagnostics_energy_start) << std::endl;
					std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-simulationSWE->diagnostics_mass_start)/simulationSWE->diagnostics_mass_start) << std::endl;
					std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-simulationSWE->diagnostics_potential_entrophy_start)/simulationSWE->diagnostics_potential_entrophy_start) << std::endl;

					if (simVars.misc.compute_errors)
					{
						std::cout << "DIAGNOSTICS BENCHMARK DIFF H:\t" << simulationSWE->benchmark.t0_error_max_abs_h_pert << std::endl;
						std::cout << "DIAGNOSTICS BENCHMARK DIFF U:\t" << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
						std::cout << "DIAGNOSTICS BENCHMARK DIFF V:\t" << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
					}

					std::cout << "[MULE] error_end_linf_h_pert: " << simulationSWE->benchmark.t0_error_max_abs_h_pert << std::endl;
					std::cout << "[MULE] error_end_linf_u: " << simulationSWE->benchmark.t0_error_max_abs_u << std::endl;
					std::cout << "[MULE] error_end_linf_v: " << simulationSWE->benchmark.t0_error_max_abs_v << std::endl;
					std::cout << std::endl;
				}



				if (simulationSWE->compute_error_to_analytical_solution)
				{
					std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << simulationSWE->benchmark.analytical_error_rms_h << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << simulationSWE->benchmark.analytical_error_rms_u << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << simulationSWE->benchmark.analytical_error_rms_v << std::endl;

					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << simulationSWE->benchmark.analytical_error_maxabs_h << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << simulationSWE->benchmark.analytical_error_maxabs_u << std::endl;
					std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << simulationSWE->benchmark.analytical_error_maxabs_v << std::endl;
				}
			}

			std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;


			delete simulationSWE;
		} // end of gui not enabled

		SimulationBenchmarkTimings::getInstance().main.stop();

		{
			std::cout << std::endl;
			SimulationBenchmarkTimings::getInstance().output();
		}
	}

	return 0;
}

