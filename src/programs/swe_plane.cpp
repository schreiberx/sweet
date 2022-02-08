/*
 * Author: Martin Schreiber <SchreiberX@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_plane_timeintegrators
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_plane_benchmarks
 */


#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif

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
#include <sweet/SWEETError.hpp>
#include <ostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#include <sweet/SimulationBenchmarkTiming.hpp>

#include "swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"
#include "swe_plane_timeintegrators/SWE_Plane_TimeSteppers.hpp"
#include "swe_plane_timeintegrators/SWE_Plane_Normal_Modes.hpp"


// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



#ifndef SWEET_MPI
#	define SWEET_MPI 1
#endif

#if SWEET_MPI
#	include <mpi.h>
#endif

#ifndef SWEET_PARAREAL
#	define SWEET_PARAREAL 1
#endif

#if SWEET_PARAREAL
	#include <parareal/Parareal.hpp>
#endif

/// general parameters
SimulationVariables simVars;

class SimulationInstance
#if SWEET_PARAREAL
		:
		public Parareal_SimulationInstance
#endif
{
public:
	// Prognostic variables
	// h: surface height (perturbation)
	// u: velocity in x-direction
	// v: velocity in y-direction
	PlaneData prog_h_pert, prog_u, prog_v;


#if SWEET_GUI
	//visualization variable
	PlaneData vis;
#endif

	// Initial values for comparison with analytical solution
	PlaneData t0_prog_h_pert, t0_prog_u, t0_prog_v;

	// Forcings
	PlaneData force_h_pert, force_u, force_v;

	PlaneDataGridMapping gridMapping;

	// Implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// Implementation of different time steppers
	SWE_Plane_TimeSteppers timeSteppersCoarse;
#endif

	//Swe benchmarks
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

	class NormalModesData
	{
	public:
		// Diagnostic information about the projection to 
		//    the linear normal wave mode eigenspace (see SWE_bench_NormalModes.hpp)
		
		PlaneData geo;    //Coefficients multiplying geostrophic mode
		PlaneData igwest; //Coefficients multiplying west gravity mode
		PlaneData igeast; //Coefficients multiplying east gravity mode
		double norm_spec;
		
	public:
		NormalModesData(
			PlaneDataConfig *planeDataConfig
		)	:
			geo(planeDataConfig),
			igwest(planeDataConfig),
			igeast(planeDataConfig)
		{
		}
	}; 

	NormalModesData normalmodes;

	// Finite difference operators
	PlaneOperators op;

	/// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;
	double diagnostics_mass_start = 0;
	double diagnostics_potential_enstrophy_start = 0;


	bool compute_error_difference_to_initial_condition = false;
	bool compute_error_to_analytical_solution = false;
	bool compute_normal_modes = false;


public:
	SimulationInstance()	:
		// Constructor to initialize the class - all variables in the SW are setup

		// Variable dimensions (mem. allocation)
		prog_h_pert(planeDataConfig),
		prog_u(planeDataConfig),
		prog_v(planeDataConfig),

#if SWEET_GUI
		vis(planeDataConfig),
#endif

		t0_prog_h_pert(planeDataConfig),
		t0_prog_u(planeDataConfig),
		t0_prog_v(planeDataConfig),

		force_h_pert(planeDataConfig),
		force_u(planeDataConfig),
		force_v(planeDataConfig),

		normalmodes(planeDataConfig),

		// Initialises operators
		op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs)
#if SWEET_PARAREAL != 0
		,
		_parareal_data_start_h(planeDataConfig), _parareal_data_start_u(planeDataConfig), _parareal_data_start_v(planeDataConfig),
		_parareal_data_fine_h(planeDataConfig), _parareal_data_fine_u(planeDataConfig), _parareal_data_fine_v(planeDataConfig),
		_parareal_data_coarse_h(planeDataConfig), _parareal_data_coarse_u(planeDataConfig), _parareal_data_coarse_v(planeDataConfig),
		_parareal_data_output_h(planeDataConfig), _parareal_data_output_u(planeDataConfig), _parareal_data_output_v(planeDataConfig),
		_parareal_data_error_h(planeDataConfig), _parareal_data_error_u(planeDataConfig), _parareal_data_error_v(planeDataConfig),
		// For parareal_SL: store penult time step in the current simulation (to be transmitted to the following time slice);
		// and penult time step received from previous time slice
		_parareal_data_coarse_previous_timestep_h(planeDataConfig), _parareal_data_coarse_previous_timestep_u(planeDataConfig), _parareal_data_coarse_previous_timestep_v(planeDataConfig),
		_parareal_data_coarse_previous_time_slice_h(planeDataConfig), _parareal_data_coarse_previous_time_slice_u(planeDataConfig), _parareal_data_coarse_previous_time_slice_v(planeDataConfig),
		// Same thing, but in the case where fine solver = SL
		_parareal_data_fine_previous_timestep_h(planeDataConfig), _parareal_data_fine_previous_timestep_u(planeDataConfig), _parareal_data_fine_previous_timestep_v(planeDataConfig),
		_parareal_data_fine_previous_time_slice_h(planeDataConfig), _parareal_data_fine_previous_time_slice_u(planeDataConfig), _parareal_data_fine_previous_time_slice_v(planeDataConfig)

#if SWEET_DEBUG
		,
		_parareal_data_fine_exact_h(planeDataConfig), _parareal_data_fine_exact_u(planeDataConfig), _parareal_data_fine_exact_v(planeDataConfig)
#endif

#endif
	{
		// Calls initialisation of the run (e.g. sets u, v, h)
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
		SimulationBenchmarkTimings::getInstance().main_setup.start();

		simVars.reset();

		if (simVars.benchmark.benchmark_name == "")
		{
			std::cout << "Benchmark scenario not selected (option --benchmark-name [string])" << std::endl;
			swePlaneBenchmarks.printBenchmarkInformation();
			SWEETError("Benchmark name not given");
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
			SWEETError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (simVars.disc.space_grid_use_c_staggering ||  !simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
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

		if (simVars.benchmark.benchmark_name == "normalmodes" )
			compute_normal_modes = true;

		update_normal_modes();
		update_diagnostics();
		
		diagnostics_energy_start = simVars.diag.total_energy;
		diagnostics_mass_start = simVars.diag.total_mass;
		diagnostics_potential_enstrophy_start = simVars.diag.total_potential_enstrophy;

		timestep_do_output();

		SimulationBenchmarkTimings::getInstance().main_setup.stop();

	}

	//Update diagnostic variables related to normal modes
	void update_normal_modes()
	{
		if (!compute_normal_modes)
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		//Setup diagnostics for normal mode projection
		SWE_bench_NormalModes::convert_allspectralmodes_to_normalmodes(
			prog_h_pert, prog_u, prog_v, simVars, // Input fields
			normalmodes.geo, normalmodes.igwest, normalmodes.igeast//Projected normal modes
		);
		
		if (simVars.timecontrol.current_timestep_nr == 0){
			//save the reference normalization parameter
			std::cout << normalmodes.geo.reduce_rms_spec() << std::endl;
			std::cout << normalmodes.igwest.reduce_rms_spec() << std::endl;
			std::cout << normalmodes.igeast.reduce_rms_spec() << std::endl;

			normalmodes.norm_spec = normalmodes.geo.reduce_sum_sq_spec()+
				normalmodes.igwest.reduce_sum_sq_spec()+
				normalmodes.igeast.reduce_sum_sq_spec();

			normalmodes.norm_spec=std::sqrt(normalmodes.norm_spec);

			if(normalmodes.norm_spec < 10e-14 ){
				normalmodes.norm_spec = 1.0;
				return;
			}
				
		}
		normalmodes.geo = normalmodes.geo / normalmodes.norm_spec;
		normalmodes.igwest = normalmodes.igwest / normalmodes.norm_spec;
		normalmodes.igeast = normalmodes.igeast / normalmodes.norm_spec;
#endif
	}


	//Update diagnostic variables related to normal modes
	void dump_normal_modes()
	{
		if (!compute_normal_modes )
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE

		SWE_bench_NormalModes::dump_all_normal_modes(simVars, normalmodes.geo, normalmodes.igwest, normalmodes.igeast);
		//normalmodes.geo.print_spectralIndex();
		//std::cout<<SWE_bench_NormalModes::bcasename <<std::endl;
#endif
		
		return;
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
			SWEETError("Max simulation time exceeded!");

#if !SWEET_PARAREAL
		timestep_do_output();
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
		char buffer[1024];

		const char* filename_template = simVars.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, simVars.timecontrol.current_simulation_time*simVars.iodata.output_time_scale);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}

	/**
	 * Write current time step info to file
	 */
	
	std::string write_output_file(
			std::stringstream &buffer
		)
	{
		const char* filename_template = "output_diag_evol.txt";
		std::ofstream file(filename_template, std::ofstream::out | std::ofstream::app);
		file << std::setprecision(12);
  		file << buffer.str() << std::endl;

		return buffer.str();
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

		if(compute_normal_modes)
			update_normal_modes();

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

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			if(compute_normal_modes){
				output_filenames += ";" + write_file_spec(normalmodes.geo, "nm_geo");
				output_filenames += ";" + write_file_spec(normalmodes.igwest, "nm_igwest");
				output_filenames += ";" + write_file_spec(normalmodes.igeast, "nm_igeast");
			}
#endif
			

		}

		if (simVars.misc.verbosity > 0)
		{
			update_diagnostics();
			compute_errors();

			if (simVars.timecontrol.current_timestep_nr == 0)
				simVars.diag.outputConfig();

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
			rows << "\t" << std::abs((simVars.diag.total_potential_enstrophy-diagnostics_potential_enstrophy_start)/diagnostics_potential_enstrophy_start);
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

#if 1
			// Normal mode stuff
			if (compute_normal_modes)
			{
				// normal modes energy
				if (simVars.timecontrol.current_timestep_nr == 0)
					header << "\tNM_GEO_RMS\tNM_IGWEST_RMS\tNM_IGEAST_RMS";

				rows << "\t" << normalmodes.geo.reduce_rms() << "\t" << normalmodes.igwest.reduce_rms_spec() << "\t" << normalmodes.igeast.reduce_rms_spec();

				//Dump to file all normal mode evolution
				dump_normal_modes();
			}
#endif

			// screen output
			if (simVars.timecontrol.current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;

#if 1
			// output to file
			if (simVars.timecontrol.current_timestep_nr == 0)
				write_output_file(header);
			write_output_file(rows);
#endif

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


#if SWEET_GUI

	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		if (simVars.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
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
	VisStuff vis_arrays[3] =
	{
			{&prog_h_pert,	"h"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
	};


	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		if (simVars.misc.vis_id < 0)
		{
			// Analytical solution at specific time on A-grid
			PlaneData ts_h_pert = t0_prog_h_pert;
			PlaneData ts_u = t0_prog_u;
			PlaneData ts_v = t0_prog_v;

			if(simVars.misc.vis_id == -1 || simVars.misc.vis_id == -2 )
			{
				//Missing to setup REXIFunctions, so invalid phi function set
				
				// Run exact solution for linear case				
				timeSteppers.l_direct->run_timestep(
						ts_h_pert, ts_u, ts_v,
						simVars.timecontrol.current_simulation_time,
						0			// initial condition given at time 0
				);
				
			}

			switch(simVars.misc.vis_id)
			{
			case -1:
				vis = ts_h_pert+simVars.sim.h0;			//Exact linear solution
				break;

			case -2:
				vis = ts_h_pert-prog_h_pert;	// difference to exact linear solution
				break;

			case -3:
				vis = t0_prog_h_pert-prog_h_pert;	// difference to initial condition
				break;

			case -4:
				vis = op.diff_c_x(prog_v) - op.diff_c_y(prog_u);	// relative vorticity
				break;
			case -5:
				vis = prog_h_pert;			//Perturbation of depth
				break;
			case -6:
				vis = normalmodes.geo ;	// geostrophic mode
				break;
			case -7:
				vis = normalmodes.igwest;	// inertia grav mode west
				break;
			case -8:
				vis = normalmodes.igeast;	// inertia grav mode east
				break;
			}

			*o_dataArray = &vis;
			*o_aspect_ratio = simVars.sim.plane_domain_size[1] / simVars.sim.plane_domain_size[0];
			return;
		}

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		if (id == 0)
		{
			vis = *vis_arrays[id].data+simVars.sim.h0;
			*o_dataArray = &vis;
		}
		else
			*o_dataArray = vis_arrays[id].data;

		vis=**o_dataArray;
		*o_aspect_ratio = simVars.sim.plane_domain_size[1] / simVars.sim.plane_domain_size[0];
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		// first, update diagnostics if required
		update_normal_modes();
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
				description = "Direct solution for h (linear only)";
				break;

			case -2:
				description = "Diff in h to exact linear solution";
				break;

			case -3:
				description = "Diff in h to initial condition";
				break;
			case -4:
				description = "Relative vorticity";
				break;
			case -5:
				description = "Depth perturbation";
				break;
			case -6:
				description = "Geostrophic wave";
				break;
			case -7:
				description = "Inertia-gravity west wave";
				break;
			case -8:
				description = "Inertia-gravity east wave";
				break;
			}
		}

		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string, "Time (days): %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
				simVars.timecontrol.current_simulation_time,
				simVars.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				simVars.timecontrol.current_timestep_nr,
				simVars.timecontrol.current_timestep_size,
				description,
				simVars.diag.total_mass,
				simVars.diag.total_energy,
				simVars.diag.total_potential_enstrophy,
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
			prog_h_pert.file_physical_saveData_ascii("swe_plane_dump_h.csv");
			prog_u.file_physical_saveData_ascii("swe_plane_dump_u.csv");
			prog_v.file_physical_saveData_ascii("swe_plane_dump_v.csv");
			break;

		case 'C':
			// dump data arrays to VTK
			prog_h_pert.file_physical_saveData_vtk("swe_plane_dump_h.vtk", "Height");
			prog_u.file_physical_saveData_vtk("swe_plane_dump_u.vtk", "U-Velocity");
			prog_v.file_physical_saveData_vtk("swe_plane_dump_v.vtk", "V-Velocity");
			break;

		case 'l':
			// load data arrays
			prog_h_pert.file_physical_loadData("swe_plane_dump_h.csv", simVars.iodata.initial_condition_input_data_binary);
			prog_u.file_physical_loadData("swe_plane_dump_u.csv", simVars.iodata.initial_condition_input_data_binary);
			prog_v.file_physical_loadData("swe_plane_dump_v.csv", simVars.iodata.initial_condition_input_data_binary);
			break;
		}
	}
#endif


	bool instability_detected()
	{
		return !(	prog_h_pert.reduce_boolean_all_finite() &&
					prog_u.reduce_boolean_all_finite() &&
					prog_v.reduce_boolean_all_finite()
				);
	}


#if SWEET_PARAREAL

	/******************************************************
	 ******************************************************
	 *       ************** PARAREAL **************
	 ******************************************************
	 ******************************************************/

	PlaneData _parareal_data_start_h, _parareal_data_start_u, _parareal_data_start_v;
	Parareal_Data_PlaneData<3> parareal_data_start;

	PlaneData _parareal_data_fine_h, _parareal_data_fine_u, _parareal_data_fine_v;
	Parareal_Data_PlaneData<3> parareal_data_fine;

	PlaneData _parareal_data_coarse_h, _parareal_data_coarse_u, _parareal_data_coarse_v;
	Parareal_Data_PlaneData<3> parareal_data_coarse;

	PlaneData _parareal_data_output_h, _parareal_data_output_u, _parareal_data_output_v;
	Parareal_Data_PlaneData<3> parareal_data_output;

	PlaneData _parareal_data_error_h, _parareal_data_error_u, _parareal_data_error_v;
	Parareal_Data_PlaneData<3> parareal_data_error;

	PlaneData _parareal_data_coarse_previous_timestep_h, _parareal_data_coarse_previous_timestep_u, _parareal_data_coarse_previous_timestep_v;
	Parareal_Data_PlaneData<3> parareal_data_coarse_previous_timestep;

	PlaneData _parareal_data_coarse_previous_time_slice_h, _parareal_data_coarse_previous_time_slice_u, _parareal_data_coarse_previous_time_slice_v;
	Parareal_Data_PlaneData<3> parareal_data_coarse_previous_time_slice;

	PlaneData _parareal_data_fine_previous_timestep_h, _parareal_data_fine_previous_timestep_u, _parareal_data_fine_previous_timestep_v;
	Parareal_Data_PlaneData<3> parareal_data_fine_previous_timestep;

	PlaneData _parareal_data_fine_previous_time_slice_h, _parareal_data_fine_previous_time_slice_u, _parareal_data_fine_previous_time_slice_v;
	Parareal_Data_PlaneData<3> parareal_data_fine_previous_time_slice;

#if SWEET_DEBUG
	PlaneData _parareal_data_fine_exact_h, _parareal_data_fine_exact_u, _parareal_data_fine_exact_v;
	Parareal_Data_PlaneData<3> parareal_data_fine_exact;
#endif

	double timeframe_start = -1;
	double timeframe_end = -1;

	bool output_data_valid = false;

	void parareal_setup()
	{
		{
			PlaneData* data_array[3] = {&_parareal_data_start_h, &_parareal_data_start_u, &_parareal_data_start_v};
			parareal_data_start.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_fine_h, &_parareal_data_fine_u, &_parareal_data_fine_v};
			parareal_data_fine.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_coarse_h, &_parareal_data_coarse_u, &_parareal_data_coarse_v};
			parareal_data_coarse.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_output_h, &_parareal_data_output_u, &_parareal_data_output_v};
			parareal_data_output.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_error_h, &_parareal_data_error_u, &_parareal_data_error_v};
			parareal_data_error.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_coarse_previous_timestep_h, &_parareal_data_coarse_previous_timestep_u, &_parareal_data_coarse_previous_timestep_v};
			parareal_data_coarse_previous_timestep.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_coarse_previous_time_slice_h, &_parareal_data_coarse_previous_time_slice_u, &_parareal_data_coarse_previous_time_slice_v};
			parareal_data_coarse_previous_time_slice.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_fine_previous_timestep_h, &_parareal_data_fine_previous_timestep_u, &_parareal_data_fine_previous_timestep_v};
			parareal_data_fine_previous_timestep.setup(data_array);
		}

		{
			PlaneData* data_array[3] = {&_parareal_data_fine_previous_time_slice_h, &_parareal_data_fine_previous_time_slice_u, &_parareal_data_fine_previous_time_slice_v};
			parareal_data_fine_previous_time_slice.setup(data_array);
		}

#if SWEET_DEBUG
		{
			PlaneData* data_array[3] = {&_parareal_data_fine_exact_h, &_parareal_data_fine_exact_u, &_parareal_data_fine_exact_v};
			parareal_data_fine_exact.setup(data_array);
		}
#endif


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
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_Data& get_reference_to_data_timestep_fine()
	{
		return parareal_data_fine;
	}

	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_Data& get_reference_to_data_timestep_coarse()
	{
		return parareal_data_coarse;
	}

	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_Data& get_reference_to_output_data()
	{
		return parareal_data_output;
	}

	/**
	 * return the penult time step of the coarse propagation
	 */
	Parareal_Data& get_reference_to_data_timestep_coarse_previous_timestep()
	{
		return parareal_data_coarse_previous_timestep;
	}


	/**
	 * return the penult time step of the fine propagation
	 */
	Parareal_Data& get_reference_to_data_timestep_fine_previous_timestep()
	{
		return parareal_data_fine_previous_timestep;
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

		*parareal_data_start.data_arrays[0] = prog_h_pert;
		*parareal_data_start.data_arrays[1] = prog_u;
		*parareal_data_start.data_arrays[2] = prog_v;

                // Useful in the first timestep of each time slice
		*parareal_data_coarse_previous_time_slice.data_arrays[0] = prog_h_pert;
		*parareal_data_coarse_previous_time_slice.data_arrays[1] = prog_u;
		*parareal_data_coarse_previous_time_slice.data_arrays[2] = prog_v;

                // Useful in the first timestep of each time slice
		*parareal_data_fine_previous_time_slice.data_arrays[0] = prog_h_pert;
		*parareal_data_fine_previous_time_slice.data_arrays[1] = prog_u;
		*parareal_data_fine_previous_time_slice.data_arrays[2] = prog_v;

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
	 * Set solution of penult coarse timestep of previous time slice
	 */
	void sim_set_data_coarse_previous_time_slice(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_coarse_previous_time_slice()" << std::endl;

		// copy to buffers
		parareal_data_coarse_previous_time_slice = i_pararealData;
	}

	/**
	 * Set solution of penult fine timestep of previous time slice
	 */
	void sim_set_data_fine_previous_time_slice(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_fine_previous_time_slice()" << std::endl;

		// copy to buffers
		parareal_data_fine_previous_time_slice = i_pararealData;
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

		prog_h_pert = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// If fine solver = SL, send penult fine time step of previous slice, except if it is the first time slice
		if ( ! simVars.disc.timestepping_method.compare("l_cn_na_sl_nd_settls") ||
                     ! simVars.disc.timestepping_method.compare("l_rexi_na_sl_nd_etdrk") ||
                     ! simVars.disc.timestepping_method.compare("l_rexi_na_sl_nd_settls") )
		{
			PlaneData h_prev = *parareal_data_fine_previous_time_slice.data_arrays[0];
			PlaneData u_prev = *parareal_data_fine_previous_time_slice.data_arrays[1];
			PlaneData v_prev = *parareal_data_fine_previous_time_slice.data_arrays[2];
			timeSteppers.master->set_previous_solution(h_prev, u_prev, v_prev);
		}

		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{
			run_timestep();
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_fine.data_arrays[0] = prog_h_pert;
		*parareal_data_fine.data_arrays[1] = prog_u;
		*parareal_data_fine.data_arrays[2] = prog_v;
	}


	/**
	 * return the data after running computations with the fine timestepping:
	 * return Y^F
	 */
	Parareal_Data& get_data_timestep_fine()
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

		prog_h_pert = *parareal_data_start.data_arrays[0];
		prog_u = *parareal_data_start.data_arrays[1];
		prog_v = *parareal_data_start.data_arrays[2];

		// reset simulation time
		simVars.timecontrol.current_simulation_time = timeframe_start;
		simVars.timecontrol.max_simulation_time = timeframe_end;
		simVars.timecontrol.current_timestep_nr = 0;

		// If coarse solver = SL, send penult coarse time step of previous slice, except if it is the first time slice
		if ( ! simVars.parareal.coarse_timestepping_method.compare("l_cn_na_sl_nd_settls") ||
                     ! simVars.parareal.coarse_timestepping_method.compare("l_rexi_na_sl_nd_etdrk") ||
                     ! simVars.parareal.coarse_timestepping_method.compare("l_rexi_na_sl_nd_settls") )
		{
			PlaneData h_prev = *parareal_data_coarse_previous_time_slice.data_arrays[0];
			PlaneData u_prev = *parareal_data_coarse_previous_time_slice.data_arrays[1];
			PlaneData v_prev = *parareal_data_coarse_previous_time_slice.data_arrays[2];
			timeSteppersCoarse.master->set_previous_solution(h_prev, u_prev, v_prev);
		}

                // Considering the case coarse timestep != time slice length
		while (simVars.timecontrol.current_simulation_time != timeframe_end)
		{

			// store previous time step
			// to be used as n-1 in SL in the next time slice
			*parareal_data_coarse_previous_timestep.data_arrays[0] = prog_h_pert;
			*parareal_data_coarse_previous_timestep.data_arrays[1] = prog_u;
			*parareal_data_coarse_previous_timestep.data_arrays[2] = prog_v;

			// allowing coarse timestepping != fine timesteppin
			timeSteppersCoarse.master->run_timestep(
								prog_h_pert, prog_u, prog_v,
								//timeframe_end - timeframe_start,
								simVars.parareal.coarse_timestep_size,
								simVars.timecontrol.current_timestep_nr++
								//simVars.timecontrol.current_timestep_nr
								//-1
			);
			simVars.timecontrol.current_simulation_time += simVars.parareal.coarse_timestep_size;
			assert(simVars.timecontrol.current_simulation_time <= timeframe_end);
		}

		// copy to buffers
		*parareal_data_coarse.data_arrays[0] = prog_h_pert;
		*parareal_data_coarse.data_arrays[1] = prog_u;
		*parareal_data_coarse.data_arrays[2] = prog_v;
	}


	/**
	 * return the solution after the coarse timestepping:
	 * return Y^C
	 */
	Parareal_Data& get_data_timestep_coarse()
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

		for (int k = 0; k < 3; k++)
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
		double convergence = -1;

		if (!i_compute_convergence_test || !output_data_valid)
		{
			for (int k = 0; k < 3; k++) {
				*parareal_data_output.data_arrays[k] = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];
				//std::cout << timeframe_end << " " << k << " " << (*parareal_data_output.data_arrays[k] - *parareal_data_fine.data_arrays[k]).reduce_maxAbs() << std::endl;
			}

			// The following lines are necessary for correctly computing the 1st iteration
			// (else, the first time step is not changed from the 0th to the 1st iteration)
			// Why?
			simVars.timecontrol.current_simulation_time = timeframe_end;
			prog_h_pert = *parareal_data_output.data_arrays[0];
			prog_u = *parareal_data_output.data_arrays[1];
			prog_v = *parareal_data_output.data_arrays[2];

			output_data_valid = true;
			return convergence;
		}

		for (int k = 0; k < 3; k++)
		{
			PlaneData tmp = *parareal_data_coarse.data_arrays[k] + *parareal_data_error.data_arrays[k];

			convergence = std::max(
					convergence,
					(*parareal_data_output.data_arrays[k]-tmp).reduce_maxAbs()
				);

			*parareal_data_output.data_arrays[k] = tmp;
		}

		simVars.timecontrol.current_simulation_time = timeframe_end;
		prog_h_pert = *parareal_data_output.data_arrays[0];
		prog_u = *parareal_data_output.data_arrays[1];
		prog_v = *parareal_data_output.data_arrays[2];

		if (compute_error_to_analytical_solution)
		{
			if (simVars.misc.compute_errors > 0)
			{
				compute_errors();
				std::cout << "maxabs error compared to analytical solution: " << benchmark.analytical_error_maxabs_h << std::endl;
			}
		}

		output_data_valid = true;
		return convergence;
	}

	/**
	 * Return the data to be forwarded to the next coarse time step interval:
	 * return Y^O
	 */
	Parareal_Data& get_output_data()
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
		Parareal_Data_PlaneData<3>& data = (Parareal_Data_PlaneData<3>&)i_data;

                // save same file but naming as slice_iter for visualizing in paraview
		std::ostringstream ss2;
		ss2 << "output_slice" << time_slice_id << "_iter" << iteration_id << ".vtk";

		std::string filename2 = ss2.str();

		data.data_arrays[0]->file_physical_saveData_vtk(filename2.c_str(), filename2.c_str());

		// save .csv files at each time step and iteration

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
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

		// Dump  data in csv, if output filename is not empty
		if (simVars.iodata.output_file_name.size() > 0)
		{
			output_filenames = "";

			output_filenames = write_file_parareal(t_h, "prog_h_pert", iteration_id);
			output_filenames += ";" + write_file_parareal(t_u, "prog_u", iteration_id);
			output_filenames += ";" + write_file_parareal(t_v, "prog_v", iteration_id);

			output_filenames += ";" + write_file_parareal(op.ke(t_u,t_v),"diag_ke", iteration_id);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			output_filenames += ";" + write_file_spec_parareal(op.ke(t_u,t_v),"diag_ke_spec", iteration_id);
#endif

			output_filenames += ";" + write_file_parareal(op.vort(t_u, t_v), "diag_vort", iteration_id);
			output_filenames += ";" + write_file_parareal(op.div(t_u, t_v), "diag_div", iteration_id);

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			if(compute_normal_modes){
				output_filenames += ";" + write_file_spec_parareal(normalmodes.geo, "nm_geo", iteration_id);
				output_filenames += ";" + write_file_spec_parareal(normalmodes.igwest, "nm_igwest", iteration_id);
				output_filenames += ";" + write_file_spec_parareal(normalmodes.igeast, "nm_igeast", iteration_id);
			}
#endif
			
		}

	}


	/**
	 * Write file to data and return string of file name (parareal)
	 */
	std::string write_file_parareal(
			const PlaneData &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, timeframe_end, iteration_id);
		i_planeData.file_physical_saveData_ascii(buffer);
		return buffer;
	}


	/**
	 * Write current time step info to file (parareal)
	 */
	
	std::string write_output_file_parareal(
			std::stringstream &buffer
		)
	{
		const char* filename_template = "output_diag_evol.txt";
		std::ofstream file(filename_template, std::ofstream::out | std::ofstream::app);
		file << std::setprecision(12);
  		file << buffer.str() << std::endl;

		return buffer.str();
	}


	/**
	 * Write spectrum info to data and return string of file name (parareal)
	 */

#if SWEET_USE_PLANE_SPECTRAL_SPACE
	std::string write_file_spec_parareal(
			const PlaneData &i_planeData,
			const char* i_name,	///< name of output variable
			int iteration_id
		)
	{
		char buffer[1024];

		const char* filename_template = "output_%s_t%020.8f_iter%03d.csv";
		sprintf(buffer, filename_template, i_name, timeframe_end, iteration_id);
		i_planeData.file_spectral_abs_saveData_ascii(buffer);
		return buffer;
	}
#endif

	void output_data_console(
			const Parareal_Data& i_data,
			int iteration_id,
			int time_slice_id
	)
	{
	}


	// check for nan in parareal (to avoid unnecessary computation)
	void check_for_nan_parareal()
	{
		int physical_size_x = parareal_data_output.data_arrays[0]->planeDataConfig->physical_data_size[0];
		int physical_size_y = parareal_data_output.data_arrays[0]->planeDataConfig->physical_data_size[1];
		for (int m = 0; m < 3; ++m)
			for (int ix = 0; ix < physical_size_x; ++ix)
				for (int iy = 0; iy < physical_size_y; ++iy)
					if ( std::isnan(parareal_data_output.data_arrays[m]->p_physical_get(ix, iy)))
						SWEETError("Instability detected in parareal!");
	}

#if SWEET_DEBUG
	/**
	* Store exact solution (full fine simulation) at the end of the time slice
	*/
	void sim_set_data_fine_exact(
			Parareal_Data &i_pararealData
	)
	{
		if (simVars.parareal.verbosity > 2)
			std::cout << "sim_set_data_fine_exact()" << std::endl;

		// copy to buffers
		parareal_data_fine_exact = i_pararealData;
	}

	/**
	* Check if solution at time k (end of time slice k-1) is exact (= fine) at iteration k
	*/
	virtual void compare_to_fine_exact()
	{
		double error = -1e10;
		double eps = 1e-10;
		for (int k = 0; k < 3; ++k)
			error = std::max(error,
					(*parareal_data_output.data_arrays[k] - *parareal_data_fine_exact.data_arrays[k]).reduce_maxAbs());

		std::cout << "Error between parareal and fine (exact) solution at t = " << timeframe_end << ": " << error << std::endl;
		if (error < eps)
			std::cout << "Parareal solution computed correctly" << std::endl;
		else
			SWEETError("Parareal solution has not been correctly computed");

	}

#endif



#endif
};



int main(int i_argc, char *i_argv[])
{
#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING_SPACE
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
			SWEETError("MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye...");
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

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


#if SWEET_PARAREAL
		simVars.parareal.printOptions();
#endif
		return -1;
	}
	
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


#if SWEET_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// only start simulation and time stepping for first rank
	if (mpi_rank == 0)
#endif
	{
#if SWEET_MPI
		std::cout << "MPI_RANK: " << mpi_rank << std::endl;
#endif

		SimulationBenchmarkTimings::getInstance().main.start();

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

#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (simVars.misc.gui_enabled)
		{
			SimulationInstance *simulationSWE = new SimulationInstance;

			SimulationBenchmarkTimings::getInstance().main_timestepping.start();

			VisSweet<SimulationInstance> visSweet(simulationSWE);

			SimulationBenchmarkTimings::getInstance().main_timestepping.stop();

			delete simulationSWE;
		}
		else
#endif
		{
			SimulationInstance *simulationSWE = new SimulationInstance;
			//Setting initial conditions and workspace - in case there is no GUI

			// also initializes diagnostics
			// already called in constructor
			//simulationSWE->reset();

#if SWEET_MPI
			MPI_Barrier(MPI_COMM_WORLD);
#endif


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
							SWEETError("INSTABILITY DETECTED");
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


#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				if (simVars.misc.verbosity > 0)
				{
					std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((simVars.diag.total_energy-simulationSWE->diagnostics_energy_start)/simulationSWE->diagnostics_energy_start) << std::endl;
					std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((simVars.diag.total_mass-simulationSWE->diagnostics_mass_start)/simulationSWE->diagnostics_mass_start) << std::endl;
					std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((simVars.diag.total_potential_enstrophy-simulationSWE->diagnostics_potential_enstrophy_start)/simulationSWE->diagnostics_potential_enstrophy_start) << std::endl;

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

#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			std::cout << std::endl;
			SimulationBenchmarkTimings::getInstance().output();
		}
	}
#if SWEET_MPI
	else	// mpi_rank != 0
	{

		if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
		{
			PlaneOperators op(planeDataConfig, simVars.sim.plane_domain_size, simVars.disc.space_use_spectral_basis_diffs);

			SWE_Plane_TS_l_rexi rexiSWE(simVars, op);

			/*
			 * Setup our little dog REXI
			 */
			rexiSWE.setup(simVars.rexi, "phi0", simVars.timecontrol.current_timestep_size);

			PlaneData prog_h_pert(planeDataConfig);
			PlaneData prog_u(planeDataConfig);
			PlaneData prog_v(planeDataConfig);

			MPI_Barrier(MPI_COMM_WORLD);

			do
			{
				// REXI time stepping
				rexiSWE.run_timestep(
						prog_h_pert,
						prog_u,
						prog_v,
						simVars.timecontrol.current_timestep_size,
						simVars.timecontrol.current_simulation_time

				);
			}
			while(!rexiSWE.final_timestep);
		}
	}

	if (simVars.disc.timestepping_method.find("rexi") != std::string::npos)
	{
		// synchronize REXI
		if (mpi_rank == 0)
			SWE_Plane_TS_l_rexi::MPI_quitWorkers(planeDataConfig);
	}

	MPI_Finalize();
#endif

	return 0;
}

