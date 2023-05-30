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
#include <sweet/core/plane/PlaneDataGridMapping.hpp>
#include "ShackPDESWEPlane_Diagnostics.hpp"
#include "benchmarks/ShackPDESWEPlaneBenchmarks.hpp"

// Benchmarks
#include "PDESWEPlane_BenchmarksCombined.hpp"

// Time steppers
#include "PDESWEPlane_TimeSteppers.hpp"

// If doing normal mode analysis
#include "PDESWEPlane_NormalModes.hpp"

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
	class SimDataAndOps
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneData_Config planeDataConfig;
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
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

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
			
			if (i_shackPlaneDataOps->space_grid_use_c_staggering)
				gridMapping.setup(i_shackPlaneDataOps, &planeDataConfig);

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
	SimDataAndOps dataAndOps;


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
		double t0_diff_error_max_abs_h_pert;
		double t0_diff_error_max_abs_u;
		double t0_diff_error_max_abs_v;

		// Error measures L2 norm
		double analytical_error_rms_h;
		double analytical_error_rms_u;
		double analytical_error_rms_v;

		// Error measures max norm
		double analytical_error_maxabs_h;
		double analytical_error_maxabs_u;
		double analytical_error_maxabs_v;

		void setup()
		{
			t0_diff_error_max_abs_h_pert = -1;
			t0_diff_error_max_abs_u = -1;
			t0_diff_error_max_abs_v = -1;

			analytical_error_rms_h = -1;
			analytical_error_rms_u = -1;
			analytical_error_rms_v = -1;

			analytical_error_maxabs_h = -1;
			analytical_error_maxabs_u = -1;
			analytical_error_maxabs_v = -1;
		}
	};

	BenchmarkErrors benchmarkErrors;

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
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWEPlaneNormalModes);
			return true;
		}

	public:
		bool setup(
			sweet::PlaneData_Config *planeDataConfig
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
		ERROR_CHECK_COND_RETURN(shackProgArgDict);
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
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		planeBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);

		pdeSWEPlaneTimeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWEPlaneTimeSteppers);

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
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_process_arguments()
	{
		shackProgArgDict.clear();
	}
	
	
	bool setup_3_main()
	{
		/*
		 * Setup Plane Data Config & Operators
		 */
		dataAndOps.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataAndOps);

		/*
		 * After we setup the plane, we can setup the time steppers and their buffers
		 */
		pdeSWEPlaneTimeSteppers.setup(&dataAndOps.ops, &shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWEPlaneTimeSteppers);

#if 0
		pdeSWEPlaneTimeSteppers.timestepper->runTimestep(
				dataAndOps.prog_h_pert, dataAndOps.prog_u, dataAndOps.prog_v,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);
#endif

#if SWEET_GUI
		vis_plane_data.setup(dataAndOps.planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		planeBenchmarksCombined.setupInitialConditions(
				dataAndOps.prog_h_pert,
				dataAndOps.prog_u,
				dataAndOps.prog_v,
				&dataAndOps.ops,
				&dataAndOps.planeDataConfig
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);

		dataAndOps.t0_prog_h_pert = dataAndOps.prog_h_pert;

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		benchmarkErrors.setup();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		shackProgArgDict.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes->setup(&dataAndOps.planeDataConfig);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*normalmodes);
		}

		if (shackPDESWEPlane->compute_errors)
		{
			//Compute difference to initial condition (makes more sense in steady state cases, but useful in others too)
			compute_error_difference_to_initial_condition = true;

			//Compute difference to analytical solution (makes more sense in linear cases, but might be useful in others too)
			compute_error_to_analytical_solution = pdeSWEPlaneTimeSteppers.linear_only;
		}
		else
		{
			compute_error_difference_to_initial_condition = false;
			compute_error_to_analytical_solution = false;
		}



		/*
		 * Load initial conditions from file if required
		 */
		if (shackIOData->initial_condition_data_filenames.size() > 0)
			dataAndOps.prog_h_pert.file_physical_loadData(shackIOData->initial_condition_data_filenames[0].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 1)
			dataAndOps.prog_u.file_physical_loadData(shackIOData->initial_condition_data_filenames[1].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 2)
			dataAndOps.prog_v.file_physical_loadData(shackIOData->initial_condition_data_filenames[2].c_str(), shackIOData->initial_condition_input_data_binary);


		if (shackPDESWEPlaneBenchmarks->benchmark_name == "normalmodes" )
			compute_normal_modes = true;

		if (compute_normal_modes)
		{
			update_normal_modes();
			update_diagnostics();
		}

		diagnostics_energy_start = shackPDESWEPlaneDiagnostics->total_energy;
		diagnostics_mass_start = shackPDESWEPlaneDiagnostics->total_mass;
		diagnostics_potential_enstrophy_start = shackPDESWEPlaneDiagnostics->total_potential_enstrophy;

		return true;
	}


	void clear_3_main()
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

		dataAndOps.clear();
	}
	

	bool setup()
	{
		if (!setup_1_registration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_main())
			return false;

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_3_main();
		clear_2_process_arguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		// keep pause state if in GUI mode
		bool run_simulation_timesteps = shackTimestepControl->run_simulation_timesteps;

		clear();

		if (!setup())
		{
			error.print();
			return false;
		}

		shackTimestepControl->run_simulation_timesteps = run_simulation_timesteps;

		return !error.exists();
	}

	// Update diagnostic variables related to normal modes
	void update_normal_modes()
	{
		if (!compute_normal_modes)
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		// Setup planeDiagnostics for normal mode projection
		normalmodes->pdeSWEPlaneNormalModes.convert_allspectralmodes_to_normalmodes(
			dataAndOps.prog_h_pert, dataAndOps.prog_u, dataAndOps.prog_v,
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
		if (dataAndOps.last_timestep_nr_update_diagnostics == shackTimestepControl->current_timestep_nr)
			return;

		dataAndOps.last_timestep_nr_update_diagnostics = shackTimestepControl->current_timestep_nr;


		if (shackPlaneDataOps->space_grid_use_c_staggering)
		{
			shackPDESWEPlaneDiagnostics->update_staggered_huv_to_mass_energy_enstrophy(
					dataAndOps.ops,
					shackPlaneDataOps,
					shackPDESWEPlane,
					dataAndOps.prog_h_pert,
					dataAndOps.prog_u,
					dataAndOps.prog_v
			);
		}
		else
		{
			shackPDESWEPlaneDiagnostics->update_nonstaggered_huv_to_mass_energy_enstrophy(
					dataAndOps.ops,
					shackPlaneDataOps,
					shackPDESWEPlane,
					dataAndOps.prog_h_pert,
					dataAndOps.prog_u,
					dataAndOps.prog_v
			);
		}
	}



	void normal_mode_analysis()
	{
		normalmodes->pdeSWEPlaneNormalModes.normal_mode_analysis(
								dataAndOps.prog_h_pert,
								dataAndOps.prog_u,
								dataAndOps.prog_v,
								3,
								&shackProgArgDict,
								this,
								&ProgramPDESWEPlane::runTimestep
						);

		std::cout << "\n Done normal mode analysis in separate class" << std::endl;
	}


	/**
	 * Execute a single simulation time step
	 */
	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();

		std::cout << "AAA " << shackTimestepControl->current_simulation_time << " " << dataAndOps.prog_h_pert.spectral_reduce_max_abs() << " " << dataAndOps.prog_h_pert.toPhys().physical_reduce_max_abs() << std::endl;

		pdeSWEPlaneTimeSteppers.timestepper->runTimestep(
				dataAndOps.prog_h_pert, dataAndOps.prog_u, dataAndOps.prog_v,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);

		shackTimestepControl->timestepHelperEnd();

#if !SWEET_PARAREAL
		timestep_do_output();
#endif
		return true;
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


		if (!shackIOData->checkDoOutput(shackTimestepControl->current_simulation_time))
			return false;

		/*
		 * File output
		 *
		 * We write everything in non-staggered output
		 */
		// For output, variables need to be on unstaggered A-grid
		sweet::PlaneData_Physical t_h(dataAndOps.planeDataConfig);
		sweet::PlaneData_Physical t_u(dataAndOps.planeDataConfig);
		sweet::PlaneData_Physical t_v(dataAndOps.planeDataConfig);

		if (shackPlaneDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
		{
			t_h = dataAndOps.prog_h_pert.toPhys();
			dataAndOps.gridMapping.mapCtoA_u(dataAndOps.prog_u.toPhys(), t_u);
			dataAndOps.gridMapping.mapCtoA_v(dataAndOps.prog_v.toPhys(), t_v);
		}
		else
		{
			t_h = dataAndOps.prog_h_pert.toPhys();
			t_u = dataAndOps.prog_u.toPhys();
			t_v = dataAndOps.prog_v.toPhys();
		}

		update_normal_modes();

		// Dump  data in csv, if output filename is not empty
		if (shackIOData->output_file_name.size() > 0)
		{
			output_filenames = "";

/////			output_filenames = write_file(t_h, "prog_h_pert");
/////			output_filenames += ";" + write_file(t_u, "prog_u");
/////			output_filenames += ";" + write_file(t_v, "prog_v");
/////
/////			output_filenames += ";" + write_file(dataAndOps.ops.ke(t_u,t_v),"diag_ke");
/////
/////#if SWEET_USE_PLANE_SPECTRAL_SPACE
/////			output_filenames += ";" + write_file_spec(dataAndOps.ops.ke(t_u,t_v),"diag_ke_spec");
/////
/////			output_filenames += ";" + write_file_spec(t_h, "prog_h_pert_spec");
/////			output_filenames += ";" + write_file_spec(t_u, "prog_u_spec");
/////			output_filenames += ";" + write_file_spec(t_v, "prog_v_spec");
/////
/////			output_filenames += ";" + write_file_spec(dataAndOps.ops.ke(t_u,t_v).toPhys(), "diag_ke_spec");
/////#endif
/////
/////			output_filenames += ";" + write_file(dataAndOps.ops.vort(t_u, t_v), "diag_vort");
/////			output_filenames += ";" + write_file(dataAndOps.ops.div(t_u, t_v), "diag_div");

			if (this->shackIOData->output_file_mode == "csv")
			{
				output_filenames = write_file(t_h, "prog_h_pert");
				output_filenames += ";" + write_file(t_u, "prog_u");
				output_filenames += ";" + write_file(t_v, "prog_v");
				output_filenames += ";" + write_file(t_h + shackPDESWEPlane->h0, "prog_h");

				output_filenames += ";" + write_file(dataAndOps.ops.vort(t_u, t_v), "diag_vrt");
				output_filenames += ";" + write_file(dataAndOps.ops.div(t_u, t_v), "diag_div");

				///output_filenames += ";" + write_file(dataAndOps.ops.ke(t_u,t_v),"diag_ke");
			}
#if SWEET_USE_PLANE_SPECTRAL_SPACE
			else
			{
				output_filenames = write_file_spec(t_h, "prog_h_pert_spec");
				output_filenames += ";" + write_file_spec(t_u, "prog_u_spec");
				output_filenames += ";" + write_file_spec(t_v, "prog_v_spec");

				output_filenames += ";" + write_file_spec(dataAndOps.ops.vort(t_u, t_v), "diag_vort_spec");
				output_filenames += ";" + write_file_spec(dataAndOps.ops.div(t_u, t_v), "diag_div_spec");

				output_filenames += ";" + write_file_spec(dataAndOps.ops.ke(t_u,t_v).toPhys(), "diag_ke_spec");
			}
#endif




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
			computeErrors();

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

				rows << "\t" << benchmarkErrors.t0_diff_error_max_abs_h_pert << "\t" << benchmarkErrors.t0_diff_error_max_abs_u << "\t" << benchmarkErrors.t0_diff_error_max_abs_v;
			}
#endif

#if 1
			if (compute_error_to_analytical_solution)
			{
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

				rows << "\t" << benchmarkErrors.analytical_error_maxabs_h << "\t" << benchmarkErrors.analytical_error_maxabs_u << "\t" << benchmarkErrors.analytical_error_maxabs_v;
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

		shackIOData->advanceNextOutput(shackTimestepControl->current_simulation_time, shackTimestepControl->max_simulation_time);
		return true;
	}


public:
	void computeErrors()
	{
		if (compute_error_difference_to_initial_condition || compute_error_to_analytical_solution)
		{
			/**
			 * Compute difference to initial condition
			 */
			if (compute_error_difference_to_initial_condition)
			{
				benchmarkErrors.t0_diff_error_max_abs_h_pert = (dataAndOps.prog_h_pert - dataAndOps.t0_prog_h_pert).toPhys().physical_reduce_max_abs();
				benchmarkErrors.t0_diff_error_max_abs_u = (dataAndOps.prog_u - dataAndOps.t0_prog_u).toPhys().physical_reduce_max_abs();
				benchmarkErrors.t0_diff_error_max_abs_v = (dataAndOps.prog_v - dataAndOps.t0_prog_v).toPhys().physical_reduce_max_abs();
			}

			// Calculate linear exact solution, if compute error requests
			if (compute_error_to_analytical_solution)
			{
				// Analytical solution at specific time on A-grid
				sweet::PlaneData_Spectral ts_h_pert = dataAndOps.t0_prog_h_pert;
				sweet::PlaneData_Spectral ts_u = dataAndOps.t0_prog_u;
				sweet::PlaneData_Spectral ts_v = dataAndOps.t0_prog_v;

				// Run exact solution for linear case
				if (pdeSWEPlaneTimeSteppers.l_direct == nullptr)
				{
					std::cerr << "Direct solution not available" << std::endl;
					exit(1);
				}

				pdeSWEPlaneTimeSteppers.l_direct->runTimestep(
						ts_h_pert, ts_u, ts_v,
						shackTimestepControl->current_simulation_time,	// time step size
						0				// initial condition given at time 0
				);

				benchmarkErrors.analytical_error_rms_h = (ts_h_pert-dataAndOps.prog_h_pert).toPhys().physical_reduce_rms();
				benchmarkErrors.analytical_error_rms_u = (ts_u-dataAndOps.prog_u).toPhys().physical_reduce_rms();
				benchmarkErrors.analytical_error_rms_v = (ts_v-dataAndOps.prog_v).toPhys().physical_reduce_rms();

				benchmarkErrors.analytical_error_maxabs_h = (ts_h_pert-dataAndOps.prog_h_pert).toPhys().physical_reduce_max_abs();
				benchmarkErrors.analytical_error_maxabs_u = (ts_u-dataAndOps.prog_u).toPhys().physical_reduce_max_abs();
				benchmarkErrors.analytical_error_maxabs_v = (ts_v-dataAndOps.prog_v).toPhys().physical_reduce_max_abs();
			}
		}
	}


public:
	void printErrors()
	{

		if (1)
		{
			update_diagnostics();

			std::cout << "DIAGNOSTICS ENERGY DIFF:\t" << std::abs((shackPDESWEPlaneDiagnostics->total_energy-diagnostics_energy_start)/diagnostics_energy_start) << std::endl;
			std::cout << "DIAGNOSTICS MASS DIFF:\t" << std::abs((shackPDESWEPlaneDiagnostics->total_mass-diagnostics_mass_start)/diagnostics_mass_start) << std::endl;
			std::cout << "DIAGNOSTICS POTENTIAL ENSTROPHY DIFF:\t" << std::abs((shackPDESWEPlaneDiagnostics->total_potential_enstrophy-diagnostics_potential_enstrophy_start)/diagnostics_potential_enstrophy_start) << std::endl;

		}

		if (shackPDESWEPlane->compute_errors)
		{
			std::cout << "BENCHMARK DIFF H(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_h_pert << std::endl;
			std::cout << "BENCHMARK DIFF U(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_u << std::endl;
			std::cout << "BENCHMARK DIFF V(t=0):\t" << benchmarkErrors.t0_diff_error_max_abs_v << std::endl;

			std::cout << "[MULE] error_end_linf_h_pert: " << benchmarkErrors.t0_diff_error_max_abs_h_pert << std::endl;
			std::cout << "[MULE] error_end_linf_u: " << benchmarkErrors.t0_diff_error_max_abs_u << std::endl;
			std::cout << "[MULE] error_end_linf_v: " << benchmarkErrors.t0_diff_error_max_abs_v << std::endl;
			std::cout << std::endl;

			std::cout << "DIAGNOSTICS ANALYTICAL RMS H:\t" << benchmarkErrors.analytical_error_rms_h << std::endl;
			std::cout << "DIAGNOSTICS ANALYTICAL RMS U:\t" << benchmarkErrors.analytical_error_rms_u << std::endl;
			std::cout << "DIAGNOSTICS ANALYTICAL RMS V:\t" << benchmarkErrors.analytical_error_rms_v << std::endl;

			std::cout << "DIAGNOSTICS ANALYTICAL MAXABS H:\t" << benchmarkErrors.analytical_error_maxabs_h << std::endl;
			std::cout << "DIAGNOSTICS ANALYTICAL MAXABS U:\t" << benchmarkErrors.analytical_error_maxabs_u << std::endl;
			std::cout << "DIAGNOSTICS ANALYTICAL MAXABS V:\t" << benchmarkErrors.analytical_error_maxabs_v << std::endl;
		}
	}


public:
	bool should_quit()
	{
		if (shackTimestepControl->max_timesteps_nr != -1)
			if (shackTimestepControl->max_timesteps_nr <= shackTimestepControl->current_timestep_nr)
				return true;

		if (shackTimestepControl->max_simulation_time != -1)
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
				runTimestep();
	}


	struct VisStuff
	{
		const sweet::PlaneData_Spectral* data;
		const char *description;
	};

	/**
	 * Arrays for online visualisation and their textual description
	 */
	VisStuff vis_arrays[3] =
	{
			{&dataAndOps.prog_h_pert,	"h"},
			{&dataAndOps.prog_u,	"u"},
			{&dataAndOps.prog_v,	"v"},
	};


	void vis_getDataArray(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		if (vis_data_id < 0)
		{
			// Analytical solution at specific time on A-grid
			sweet::PlaneData_Spectral ts_h_pert = dataAndOps.t0_prog_h_pert;
			sweet::PlaneData_Spectral ts_u = dataAndOps.t0_prog_u;
			sweet::PlaneData_Spectral ts_v = dataAndOps.t0_prog_v;

			if (vis_data_id == -1 || vis_data_id == -2 )
			{
				// Missing to setup REXIFunctions, so invalid phi function set

				if (pdeSWEPlaneTimeSteppers.l_direct)
				{
					// Run exact solution for linear case
					pdeSWEPlaneTimeSteppers.l_direct->runTimestep(
							ts_h_pert, ts_u, ts_v,
							shackTimestepControl->current_simulation_time,
							0			// initial condition given at time 0
					);
				}
				else
				{
					SWEETError("Direct solution not available");
				}
			}

			sweet::PlaneData_Spectral vis_tmp;
			switch(vis_data_id)
			{
			case -1:
				vis_tmp = ts_h_pert+shackPDESWEPlane->h0;			// Exact linear solution
				break;

			case -2:
				vis_tmp = ts_h_pert-dataAndOps.prog_h_pert;			// difference to exact linear solution
				break;

			case -3:
				vis_tmp = dataAndOps.t0_prog_h_pert-dataAndOps.prog_h_pert;	// difference to initial condition
				break;

			case -4:
				vis_tmp = dataAndOps.ops.diff_c_x(dataAndOps.prog_v) - dataAndOps.ops.diff_c_y(dataAndOps.prog_u);	// relative vorticity
				break;
			case -5:
				vis_tmp = dataAndOps.prog_h_pert;			//Perturbation of depth
				break;
			case -6:
				vis_tmp = normalmodes->geo ;	// geostrophic mode
				break;
			case -7:
				vis_tmp = normalmodes->igwest;	// inertia grav mode west
				break;
			case -8:
				vis_tmp = normalmodes->igeast;	// inertia grav mode east
				break;
			}

			vis_plane_data = vis_tmp.toPhys();
			*o_dataArray = &vis_plane_data;
			*o_aspect_ratio = shackPlaneDataOps->plane_domain_size[1] / shackPlaneDataOps->plane_domain_size[0];
			return;
		}

		int id = vis_data_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));

		sweet::PlaneData_Spectral vis_tmp;
		if (id == 0)
		{
			vis_tmp = *vis_arrays[id].data+shackPDESWEPlane->h0;
		}
		else
		{
			vis_tmp = *(vis_arrays[id].data);
		}

		vis_plane_data = vis_tmp.toPhys();
		*o_dataArray = &vis_plane_data;

		*o_aspect_ratio = shackPlaneDataOps->plane_domain_size[1] / shackPlaneDataOps->plane_domain_size[0];
	}



	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		// first, update planeDiagnostics if required
		update_normal_modes();
		update_diagnostics();

		const char* description = "";
		if (vis_data_id >= 0)
		{
			int id = vis_data_id % (sizeof(vis_arrays)/sizeof(*vis_arrays));
			description = vis_arrays[id].description;
		}
		else
		{
			switch (vis_data_id)
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
				vis_plane_data.physical_reduce_max(),
				vis_plane_data.physical_reduce_min()
			);

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
			vis_data_id++;
			break;

		case 'V':
			vis_data_id--;
			break;

#if 0
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
#endif
		}
	}
#endif


	bool instability_detected()
	{
		return !(	dataAndOps.prog_h_pert.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps.prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps.prog_v.toPhys().physical_reduce_boolean_all_finite()
				);
	}

};


#endif
