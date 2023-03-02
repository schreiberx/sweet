/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_SWEPLANE_PROGRAMPDESWEPLANE_HPP_
#define SRC_PROGRAMS_PDE_SWEPLANE_PROGRAMPDESWEPLANE_HPP_



// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/plane/Plane.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include "ShackPDESWEPlaneDiagnostics.hpp"
#include "benchmarks/ShackPDESWEPlaneBenchmarks.hpp"

// Benchmarks
#include "PDESWEPlaneBenchmarksCombined.hpp"

// Time steppers
#include "PDESWEPlaneTimeSteppers.hpp"

// If doing normal mode analysis
#include "PDESWEPlaneNormalModes.hpp"

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
#endif


#if SWEET_MPI
#	include <mpi.h>
#endif

#if SWEET_PARAREAL
	#include <sweet/parareal/Parareal.hpp>
#endif

#if SWEET_XBRAID
	#include <sweet/xbraid/XBraid_sweet_lib.hpp>
#endif


class ProgramPDESWEPlane
#if SWEET_GUI
		:	public SimulationGUICallbacks
#endif
{
public:
	sweet::ErrorBase error;

	/*
	 * Just a class to store simulation data all together
	 */
	class Data
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneDataConfig planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h_pert;
		sweet::PlaneData_Spectral prog_u;
		sweet::PlaneData_Spectral prog_v;

		// TODO: Get rid of me right here
		// Initial values for comparison with analytical solution
		sweet::PlaneData_Spectral t0_prog_h_pert;
		sweet::PlaneData_Spectral t0_prog_u;
		sweet::PlaneData_Spectral t0_prog_v;

		// Forcings
		sweet::PlaneData_Spectral force_prog_h_pert;
		sweet::PlaneData_Spectral force_prog_u;
		sweet::PlaneData_Spectral force_prog_v;

		// Mapping between grids
		sweet::PlaneDataGridMapping gridMapping;
		
		// Diagnostics measures
		int last_timestep_nr_update_diagnostics = -1;
			
		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			planeDataConfig.setupAuto(*i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(ops);

			prog_h_pert.setup(planeDataConfig);
			prog_u.setup(planeDataConfig);
			prog_v.setup(planeDataConfig);

			t0_prog_h_pert.setup(planeDataConfig);
			t0_prog_u.setup(planeDataConfig);
			t0_prog_v.setup(planeDataConfig);

			force_prog_h_pert.setup(planeDataConfig);
			force_prog_u.setup(planeDataConfig);
			force_prog_v.setup(planeDataConfig);

			last_timestep_nr_update_diagnostics = -1;
			
			return true;
		}

		void clear()
		{
			prog_h_pert.clear();
			prog_u.clear();
			prog_v.clear();
			
			t0_prog_h_pert.clear();
			t0_prog_u.clear();
			t0_prog_v.clear();
			
			force_prog_h_pert.clear();
			force_prog_u.clear();
			force_prog_v.clear();

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Simulation data
	Data data;


	// time integrators
	PDESWEPlaneTimeSteppers pdeSWEPlaneTimeSteppers;

#if SWEET_PARAREAL
	// Implementation of different time steppers
	PDESWEPlaneTimeSteppers timeSteppersCoarse;
#endif

	// Handler to all benchmarks
	PDESWEPlaneBenchmarksCombined planeBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWEPlane *shackPDESWEPlane;
	ShackPDESWEPlaneTimeDiscretization *shackTimeDisc;
	ShackPDESWEPlaneBenchmarks *shackPDESWEPlaneBenchmarks;
	ShackPDESWEPlaneDiagnostics *shackPDESWEPlaneDiagnostics;

#if SWEET_GUI
	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical vis_plane_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 0;

	// Which primitive to use for rendering
	int vis_data_id = 0;

#endif

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
		sweet::ErrorBase error;

		// Diagnostic information about the projection to
		//    the linear normal wave mode eigenspace (see SWE_bench_normalmodes->hpp)

		sweet::PlaneData_Spectral geo;    //Coefficients multiplying geostrophic mode
		sweet::PlaneData_Spectral igwest; //Coefficients multiplying west gravity mode
		sweet::PlaneData_Spectral igeast; //Coefficients multiplying east gravity mode
		double norm_spec;

		PDESWEPlaneNormalModes pdeSWEPlaneNormalModes;

	public:
		bool shackRegistration(sweet::ShackDictionary &io_dict)
		{
			pdeSWEPlaneNormalModes.shackRegistration(io_dict);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(pdeSWEPlaneNormalModes);
			return true;
		}

	public:
		bool setup(
			sweet::PlaneDataConfig *planeDataConfig
		)
		{
			geo.setup(planeDataConfig);
			igwest.setup(planeDataConfig);
			igeast.setup(planeDataConfig);

			return true;
		}

	public:
		bool clear()
		{
			geo.clear();
			igwest.clear();
			igeast.clear();
			return true;
		}
	};

	NormalModesData *normalmodes;

	/// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;
	double diagnostics_mass_start = 0;
	double diagnostics_potential_enstrophy_start = 0;


	bool compute_error_difference_to_initial_condition = false;
	bool compute_error_to_analytical_solution = false;
	bool compute_normal_modes = false;


public:
	ProgramPDESWEPlane(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackPDESWEPlane(nullptr),
		shackTimeDisc(nullptr),
		shackPDESWEPlaneBenchmarks(nullptr),
		shackPDESWEPlaneDiagnostics(nullptr)
	{
		ERROR_CHECK_WITH_RETURN(shackProgArgDict);
	}


	bool setup_1_registration()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackPDESWEPlane = shackProgArgDict.getAutoRegistration<ShackPDESWEPlane>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneTimeDiscretization>();
		shackPDESWEPlaneBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneBenchmarks>();
		shackPDESWEPlaneDiagnostics = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneDiagnostics>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		planeBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(planeBenchmarksCombined);

		pdeSWEPlaneTimeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(pdeSWEPlaneTimeSteppers);

		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes = new NormalModesData;
			normalmodes->shackRegistration(shackProgArgDict);
		}

		return true;
	}

	void clear_1_shackRegistration()
	{
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			delete normalmodes;
			normalmodes = nullptr;
		}

		shackPlaneDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		planeBenchmarksCombined.clear();
		pdeSWEPlaneTimeSteppers.clear();
		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		shackProgArgDict.setup();

		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_process_arguments()
	{
		shackProgArgDict.clear();
	}
	
	
	bool setup_3_data()
	{
		/*
		 * Setup Plane Data Config & Operators
		 */
		data.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(data);

		/*
		 * After we setup the plane, we can setup the time steppers and their buffers
		 */
		pdeSWEPlaneTimeSteppers.setup(&data.ops, &shackProgArgDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(pdeSWEPlaneTimeSteppers);

		pdeSWEPlaneTimeSteppers.master->run_timestep(
				data.prog_h_pert, data.prog_u, data.prog_v,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);

#if SWEET_GUI
		vis_plane_data.setup(data.planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		planeBenchmarksCombined.setupInitialConditions(
				data.prog_h_pert,
				data.prog_u,
				data.prog_v,
				&data.ops,
				&data.planeDataConfig
			);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(planeBenchmarksCombined);

		data.t0_prog_h_pert = data.prog_h_pert;

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		shackProgArgDict.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes->setup(&data.planeDataConfig);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(*normalmodes);
		}

		return true;
	}
	void clear_3_data()
	{
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes->clear();
			delete normalmodes;
			normalmodes = nullptr;
		}
		
#if SWEET_GUI
		vis_plane_data.clear();
#endif

		pdeSWEPlaneTimeSteppers.clear();

		data.clear();
	}
	

	bool setup()
	{
		if (!setup_1_registration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
			return false;

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_3_data();
		clear_2_process_arguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		return !error.exists();
	}
#if 0
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


		// Initialise planeDiagnostics
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

		shackTimestepControl->current_timestep_nr = 0;
		shackTimestepControl->current_simulation_time = 0;

		// set to some values for first touch NUMA policy (HPC stuff)
#if SWEET_USE_PLANE_SPECTRAL_SPACE
		prog_h_pert.spectral_set_value(std::complex<double>(0, 0));
		prog_u.spectral_set_value(std::complex<double>(0, 0));
		prog_v.spectral_set_value(std::complex<double>(0, 0));
#endif

///		// Setup prog vars
///		prog_h_pert.physical_set_value(shackPDESWEPlane->h0);
///		prog_u.physical_set_value(0);
///		prog_v.physical_set_value(0);

		// Check if input parameters are adequate for this simulation
		if (shackPlaneDataOps->space_grid_use_c_staggering && simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Staggering and spectral basis not supported!");

#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		if (shackPlaneDataOps->space_grid_use_c_staggering ||  !simVars.disc.space_use_spectral_basis_diffs)
			SWEETError("Finite differences and spectral dealisiang should not be used together! Please compile without dealiasing.");
#endif


		if (shackPlaneDataOps->space_grid_use_c_staggering)
			gridMapping.setup(simVars, planeDataConfig);

		swePlaneBenchmarks.setupInitialConditions(t0_prog_h_pert, t0_prog_u, t0_prog_v, simVars, ops);

		prog_h_pert = t0_prog_h_pert;
		prog_u = t0_prog_u;
		prog_v = t0_prog_v;

		// Load data, if requested
		if (shackIOData->initial_condition_data_filenames.size() > 0)
			prog_h_pert.file_physical_loadData(shackIOData->initial_condition_data_filenames[0].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 1)
			prog_u.file_physical_loadData(shackIOData->initial_condition_data_filenames[1].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 2)
			prog_v.file_physical_loadData(shackIOData->initial_condition_data_filenames[2].c_str(), shackIOData->initial_condition_input_data_binary);


		pdeSWEPlaneTimeSteppers.setup(
				simVars.disc.timestepping_method,
				simVars.disc.timestepping_order,
				simVars.disc.timestepping_order2,
				ops,
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
			compute_error_to_analytical_solution = pdeSWEPlaneTimeSteppers.linear_only;
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

		diagnostics_energy_start = shackPDESWEPlaneDiagnostics->total_energy;
		diagnostics_mass_start = shackPDESWEPlaneDiagnostics->total_mass;
		diagnostics_potential_enstrophy_start = shackPDESWEPlaneDiagnostics->total_potential_enstrophy;

		timestep_do_output();

		SimulationBenchmarkTimings::getInstance().main_setup.stop();

	}
#endif
	// Update diagnostic variables related to normal modes
	void update_normal_modes()
	{
		if (!compute_normal_modes)
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		// Setup planeDiagnostics for normal mode projection
		normalmodes->pdeSWEPlaneNormalModes.convert_allspectralmodes_to_normalmodes(
			data.prog_h_pert, data.prog_u, data.prog_v,
			normalmodes->geo, normalmodes->igwest, normalmodes->igeast
		);

		if (shackTimestepControl->current_timestep_nr == 0){
			//save the reference normalization parameter
			std::cout << normalmodes->geo.spectral_reduce_rms() << std::endl;
			std::cout << normalmodes->igwest.spectral_reduce_rms() << std::endl;
			std::cout << normalmodes->igeast.spectral_reduce_rms() << std::endl;

			normalmodes->norm_spec = normalmodes->geo.spectral_reduce_sum_sqr_quad()+
				normalmodes->igwest.spectral_reduce_sum_sqr_quad()+
				normalmodes->igeast.spectral_reduce_sum_sqr_quad();

			normalmodes->norm_spec=std::sqrt(normalmodes->norm_spec);

			if (normalmodes->norm_spec < 10e-14 ){
				normalmodes->norm_spec = 1.0;
				return;
			}

		}
		normalmodes->geo = normalmodes->geo / normalmodes->norm_spec;
		normalmodes->igwest = normalmodes->igwest / normalmodes->norm_spec;
		normalmodes->igeast = normalmodes->igeast / normalmodes->norm_spec;
#endif
	}


	//Update diagnostic variables related to normal modes
	void dump_normal_modes()
	{
		if (!compute_normal_modes )
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PDESWEPlaneBench_NormalModes n;
		n.dump_all_normal_modes(normalmodes->geo, normalmodes->igwest, normalmodes->igeast);
#endif

		return;
	}

	//Calculate the model diagnostics
	void update_diagnostics()
	{
		// assure, that the planeDiagnostics are only updated for new time steps
		if (data.last_timestep_nr_update_diagnostics == shackTimestepControl->current_timestep_nr)
			return;

		data.last_timestep_nr_update_diagnostics = shackTimestepControl->current_timestep_nr;


		if (shackPlaneDataOps->space_grid_use_c_staggering)
		{
			shackPDESWEPlaneDiagnostics->update_staggered_huv_to_mass_energy_enstrophy(
					data.ops,
					shackPlaneDataOps,
					shackPDESWEPlane,
					data.prog_h_pert,
					data.prog_u,
					data.prog_v
			);
		}
		else
		{
			shackPDESWEPlaneDiagnostics->update_nonstaggered_huv_to_mass_energy_enstrophy(
					data.ops,
					shackPlaneDataOps,
					shackPDESWEPlane,
					data.prog_h_pert,
					data.prog_u,
					data.prog_v
			);
		}
	}



	void normal_mode_analysis()
	{
		normalmodes->pdeSWEPlaneNormalModes.normal_mode_analysis(
								data.prog_h_pert,
								data.prog_u,
								data.prog_v,
								3,
								&shackProgArgDict,
								this,
								&ProgramPDESWEPlane::run_timestep
						);

		std::cout << "\n Done normal mode analysis in separate class" << std::endl;
	}


	/**
	 * Execute a single simulation time step
	 */
	void run_timestep()
	{
		if (shackTimestepControl->current_simulation_time + shackTimestepControl->current_timestep_size > shackTimestepControl->max_simulation_time)
			shackTimestepControl->current_timestep_size = shackTimestepControl->max_simulation_time - shackTimestepControl->current_simulation_time;

		pdeSWEPlaneTimeSteppers.master->run_timestep(
				data.prog_h_pert, data.prog_u, data.prog_v,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);


		// advance time step and provide information to parameters
		shackTimestepControl->current_simulation_time += shackTimestepControl->current_timestep_size;
		shackTimestepControl->current_timestep_nr++;

		if (shackTimestepControl->current_simulation_time > shackTimestepControl->max_simulation_time)
			SWEETError("Max simulation time exceeded!");

#if !SWEET_PARAREAL
		timestep_do_output();
#endif
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file(
			const sweet::PlaneData_Spectral &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		char buffer[1024];

		// TODO: convert spectral datato physical

		const char* filename_template = shackIOData->output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale);
		i_planeData.toPhys().file_physical_saveData_ascii(buffer);
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
			const sweet::PlaneData_Spectral &i_planeData,
			const char* i_name	///< name of output variable
		)
	{
		char buffer[1024];

		const char* filename_template = shackIOData->output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale);
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
		if (shackPDESWEPlane->normal_mode_analysis_generation > 0)
			return false;

		// output each time step
		if (shackIOData->output_each_sim_seconds < 0)
			return false;

		if (shackIOData->output_next_sim_seconds-shackIOData->output_next_sim_seconds*(1e-12) > shackTimestepControl->current_simulation_time)
			return false;

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::PlaneData_Physical t_h(data.planeDataConfig);
		sweet::PlaneData_Physical t_u(data.planeDataConfig);
		sweet::PlaneData_Physical t_v(data.planeDataConfig);

		if (shackPlaneDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			t_h = data.prog_h_pert.toPhys();
			data.gridMapping.mapCtoA_u(data.prog_u.toPhys(), t_u);
			data.gridMapping.mapCtoA_v(data.prog_v.toPhys(), t_v);
		}
		else
		{
			t_h = data.prog_h_pert.toPhys();
			t_u = data.prog_u.toPhys();
			t_v = data.prog_v.toPhys();
		}

		//std::cout << simVars.inputoutput.output_next_sim_seconds << "\t" << shackTimestepControl->current_simulation_time << std::endl;

		if (compute_normal_modes)
			update_normal_modes();

		// Dump  data in csv, if output filename is not empty
		if (shackIOData->output_file_name.size() > 0)
		{
			output_filenames = "";

			output_filenames = write_file(t_h, "prog_h_pert");
			output_filenames += ";" + write_file(t_u, "prog_u");
			output_filenames += ";" + write_file(t_v, "prog_v");

			output_filenames += ";" + write_file(data.ops.ke(t_u,t_v),"diag_ke");

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			output_filenames += ";" + write_file_spec(data.ops.ke(t_u,t_v),"diag_ke_spec");

			output_filenames += ";" + write_file_spec(t_h, "prog_h_pert_spec");
			output_filenames += ";" + write_file_spec(t_u, "prog_u_spec");
			output_filenames += ";" + write_file_spec(t_v, "prog_v_spec");

			output_filenames += ";" + write_file_spec(data.ops.ke(t_u,t_v).toPhys(), "diag_ke_spec");
#endif

			output_filenames += ";" + write_file(data.ops.vort(t_u, t_v), "diag_vort");
			output_filenames += ";" + write_file(data.ops.div(t_u, t_v), "diag_div");

#if SWEET_USE_PLANE_SPECTRAL_SPACE
			if (compute_normal_modes){
				output_filenames += ";" + write_file_spec(normalmodes->geo, "nm_geo");
				output_filenames += ";" + write_file_spec(normalmodes->igwest, "nm_igwest");
				output_filenames += ";" + write_file_spec(normalmodes->igeast, "nm_igeast");
			}
#endif
		}

		if (shackIOData->verbosity > 0)
		{
			update_diagnostics();
			compute_errors();

			std::stringstream header;
			std::stringstream rows;

			rows << std::setprecision(16);

			// Prefix
			if (shackTimestepControl->current_timestep_nr == 0)
				header << "DATA";
			rows << "DATA";

			// Time
			if (shackTimestepControl->current_timestep_nr == 0)
				header << "\tT";
			rows << "\t" << shackTimestepControl->current_simulation_time;

#if 1
			// Mass, Energy, Enstrophy
			header << "\tTOTAL_MASS\tTOTAL_ENERGY\tPOT_ENSTROPHY";
			rows << "\t" << shackPDESWEPlaneDiagnostics->total_mass;
			rows << "\t" << shackPDESWEPlaneDiagnostics->total_energy;
			rows << "\t" << shackPDESWEPlaneDiagnostics->total_potential_enstrophy;

			// Mass, Energy, Enstrophy
			header << "\tTOTAL_MASS_REL_ERROR\tTOTAL_ENERGY_REL_ERROR\tPOT_ENSTROPHY_REL_ERROR";
			rows << "\t" << std::abs((shackPDESWEPlaneDiagnostics->total_mass-diagnostics_mass_start)/diagnostics_mass_start);
			rows << "\t" << std::abs((shackPDESWEPlaneDiagnostics->total_energy-diagnostics_energy_start)/diagnostics_energy_start) ;
			rows << "\t" << std::abs((shackPDESWEPlaneDiagnostics->total_potential_enstrophy-diagnostics_potential_enstrophy_start)/diagnostics_potential_enstrophy_start);
#endif

#if 1
			if (compute_error_difference_to_initial_condition)
			{
				// Difference to initial condition
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tDIFF_MAXABS_H0\tDIFF_MAXABS_U0\tDIFF_MAXABS_V0";

				rows << "\t" << benchmark.t0_error_max_abs_h_pert << "\t" << benchmark.t0_error_max_abs_u << "\t" << benchmark.t0_error_max_abs_v;
			}
#endif

#if 1
			if (compute_error_to_analytical_solution)
			{
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

				rows << "\t" << benchmark.analytical_error_maxabs_h << "\t" << benchmark.analytical_error_maxabs_u << "\t" << benchmark.analytical_error_maxabs_v;
			}
#endif

#if 1
			// Normal mode stuff
			if (compute_normal_modes)
			{
				// normal modes energy
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tNM_GEO_RMS\tNM_IGWEST_RMS\tNM_IGEAST_RMS";

				rows << "\t" << normalmodes->geo.spectral_reduce_rms();
				rows << "\t" << normalmodes->igwest.spectral_reduce_rms();
				rows << "\t" << normalmodes->igeast.spectral_reduce_rms();

				//Dump to file all normal mode evolution
				dump_normal_modes();
			}
#endif

			// screen output
			if (shackTimestepControl->current_timestep_nr == 0)
				o_ostream << header.str() << std::endl;

			o_ostream << rows.str() << std::endl;

#if 1
			// output to file
			if (shackTimestepControl->current_timestep_nr == 0)
				write_output_file(header);
			write_output_file(rows);
#endif

#if 1
			if (diagnostics_mass_start > 0.00001 && std::abs((shackPDESWEPlaneDiagnostics->total_mass-diagnostics_mass_start)/diagnostics_mass_start) > 10000000.0)
			{
				std::cerr << "\n DIAGNOSTICS MASS DIFF TOO LARGE:\t" << std::abs((shackPDESWEPlaneDiagnostics->total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			}
#endif

		}

		if (shackIOData->output_next_sim_seconds == shackTimestepControl->max_simulation_time)
		{
			shackIOData->output_next_sim_seconds = std::numeric_limits<double>::infinity();
		}
		else
		{
			while (shackIOData->output_next_sim_seconds-shackIOData->output_next_sim_seconds*(1e-12) <= shackTimestepControl->current_simulation_time)
				shackIOData->output_next_sim_seconds += shackIOData->output_each_sim_seconds;

			if (shackIOData->output_next_sim_seconds > shackTimestepControl->max_simulation_time)
				shackIOData->output_next_sim_seconds = shackTimestepControl->max_simulation_time;
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
				benchmark.t0_error_max_abs_h_pert = (data.prog_h_pert -data. t0_prog_h_pert).toPhys().physical_reduce_max_abs();
				benchmark.t0_error_max_abs_u = (data.prog_u - data.t0_prog_u).toPhys().physical_reduce_max_abs();
				benchmark.t0_error_max_abs_v = (data.prog_v - data.t0_prog_v).toPhys().physical_reduce_max_abs();
			}

			// Calculate linear exact solution, if compute error requests
			if (compute_error_to_analytical_solution)
			{
				// Analytical solution at specific time on A-grid
				sweet::PlaneData_Spectral ts_h_pert = data.t0_prog_h_pert;
				sweet::PlaneData_Spectral ts_u = data.t0_prog_u;
				sweet::PlaneData_Spectral ts_v = data.t0_prog_v;

				// Run exact solution for linear case
				if (pdeSWEPlaneTimeSteppers.l_direct == nullptr)
				{
					std::cerr << "Direct solution not available" << std::endl;
					exit(1);
				}
				pdeSWEPlaneTimeSteppers.l_direct->run_timestep(
						ts_h_pert, ts_u, ts_v,
						shackTimestepControl->current_simulation_time,	// time step size
						0				// initial condition given at time 0
				);

				benchmark.analytical_error_rms_h = (ts_h_pert-data.prog_h_pert).toPhys().physical_reduce_rms();
				benchmark.analytical_error_rms_u = (ts_u-data.prog_u).toPhys().physical_reduce_rms();
				benchmark.analytical_error_rms_v = (ts_v-data.prog_v).toPhys().physical_reduce_rms();

				benchmark.analytical_error_maxabs_h = (ts_h_pert-data.prog_h_pert).toPhys().physical_reduce_max_abs();
				benchmark.analytical_error_maxabs_u = (ts_u-data.prog_u).toPhys().physical_reduce_max_abs();
				benchmark.analytical_error_maxabs_v = (ts_v-data.prog_v).toPhys().physical_reduce_max_abs();
			}
		}
	}



public:
	bool should_quit()
	{
		if (
				shackTimestepControl->max_timesteps_nr != -1 &&
				shackTimestepControl->max_timesteps_nr <= shackTimestepControl->current_timestep_nr
		)
			return true;

		if (!std::isinf(shackTimestepControl->max_simulation_time))
			if (shackTimestepControl->max_simulation_time <= shackTimestepControl->current_simulation_time+shackTimestepControl->max_simulation_time*1e-10)	// care about roundoff errors with 1e-10
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
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
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
	VisStuff vis_arrays[3] =
	{
			{&prog_h_pert,	"h"},
			{&prog_u,	"u"},
			{&prog_v,	"v"},
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
		if (simVars.misc.vis_id < 0)
		{
			// Analytical solution at specific time on A-grid
			sweet::PlaneData_Physical ts_h_pert = t0_prog_h_pert;
			sweet::PlaneData_Physical ts_u = t0_prog_u;
			sweet::PlaneData_Physical ts_v = t0_prog_v;

			if (simVars.misc.vis_id == -1 || simVars.misc.vis_id == -2 )
			{
				//Missing to setup REXIFunctions, so invalid phi function set

				// Run exact solution for linear case
				pdeSWEPlaneTimeSteppers.l_direct->run_timestep(
						ts_h_pert, ts_u, ts_v,
						shackTimestepControl->current_simulation_time,
						0			// initial condition given at time 0
				);

			}

			switch(simVars.misc.vis_id)
			{
			case -1:
				vis = ts_h_pert+shackPDESWEPlane->h0;			//Exact linear solution
				break;

			case -2:
				vis = ts_h_pert-prog_h_pert;	// difference to exact linear solution
				break;

			case -3:
				vis = t0_prog_h_pert-prog_h_pert;	// difference to initial condition
				break;

			case -4:
				vis = ops.diff_c_x(prog_v) - ops.diff_c_y(prog_u);	// relative vorticity
				break;
			case -5:
				vis = prog_h_pert;			//Perturbation of depth
				break;
			case -6:
				vis = normalmodes->geo ;	// geostrophic mode
				break;
			case -7:
				vis = normalmodes->igwest;	// inertia grav mode west
				break;
			case -8:
				vis = normalmodes->igeast;	// inertia grav mode east
				break;
			}

			*o_dataArray = &vis;
			*o_aspect_ratio = shackPlaneDataOps->plane_domain_size[1] / shackPlaneDataOps->plane_domain_size[0];
			return;
		}

		int id = simVars.misc.vis_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		if (id == 0)
		{
			vis = *vis_arrays[id].data+shackPDESWEPlane->h0;
			*o_dataArray = &vis;
		}
		else
			*o_dataArray = vis_arrays[id].data;

		vis=**o_dataArray;
		*o_aspect_ratio = shackPlaneDataOps->plane_domain_size[1] / shackPlaneDataOps->plane_domain_size[0];
	}



	/**
	 * return status string for window title
	 */
	const char* vis_getStatusString()
	{
		// first, update planeDiagnostics if required
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
				shackTimestepControl->current_simulation_time,
				shackTimestepControl->current_simulation_time/(60.0*60.0*24.0),
				shackTimestepControl->current_timestep_nr,
				shackTimestepControl->current_timestep_size,
				description,
				shackPDESWEPlaneDiagnostics->total_mass,
				shackPDESWEPlaneDiagnostics->total_energy,
				shackPDESWEPlaneDiagnostics->total_potential_enstrophy,
				vis.physical_reduce_max(),
				vis.physical_reduce_min() );

		return title_string;
	}



	void vis_pause()
	{
		shackTimestepControl->run_simulation_timesteps = !shackTimestepControl->run_simulation_timesteps;
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
			prog_h_pert.file_physical_loadData("swe_plane_dump_h.csv", shackIOData->initial_condition_input_data_binary);
			prog_u.file_physical_loadData("swe_plane_dump_u.csv", shackIOData->initial_condition_input_data_binary);
			prog_v.file_physical_loadData("swe_plane_dump_v.csv", shackIOData->initial_condition_input_data_binary);
			break;
		}
	}
#endif


	bool instability_detected()
	{
		return !(	data.prog_h_pert.toPhys().physical_reduce_boolean_all_finite() &&
					data.prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					data.prog_v.toPhys().physical_reduce_boolean_all_finite()
				);
	}

};


#endif
