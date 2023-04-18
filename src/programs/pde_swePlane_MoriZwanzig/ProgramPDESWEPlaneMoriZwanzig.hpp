/*
 * Author: Joao STEINSTRAESSER <joao.steinstraesser@usp.br>
 */

#ifndef SRC_PROGRAMS_PDE_SWEPLANE_MORI_ZWANZIG_PROGRAMPDESWEPLANEMORIZWANZIG_HPP_
#define SRC_PROGRAMS_PDE_SWEPLANE_MORI_ZWANZIG_PROGRAMPDESWEPLANEMORIZWANZIG_HPP_


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
#include "../pde_swePlane/ShackPDESWEPlane_Diagnostics.hpp"
///#include "../pde_swePlane/benchmarks/ShackPDESWEPlaneBenchmarks.hpp"
#include "benchmarks/ShackPDESWEPlaneMoriZwanzigBenchmarks.hpp"

// Benchmarks
////#include "../pde_swePlane/PDESWEPlane_BenchmarksCombined.hpp"
#include "PDESWEPlaneMoriZwanzig_BenchmarksCombined.hpp"

// Time steppers
#include "PDESWEPlaneMoriZwanzig_TimeSteppers.hpp"

//////// If doing normal mode analysis
//////#include "../pde_swePlane/PDESWEPlane_NormalModes.hpp"

#if SWEET_MPI
#	include <mpi.h>
#endif

#if SWEET_PARAREAL
	#include <sweet/parareal/Parareal.hpp>
#endif

#if SWEET_XBRAID
	#include <sweet/xbraid/XBraid_sweet_lib.hpp>
#endif

#include <sweet/core/StopwatchBox.hpp>

class ProgramPDESWEPlaneMoriZwanzig
{
public:
	sweet::ErrorBase error;

	PDESWEPlaneMoriZwanzigProjection projection;

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

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Simulation data
	SimDataAndOps dataAndOps_SP;
	SimDataAndOps dataAndOps_SQ;
	SimDataAndOps dataAndOps_FQ;


	// time integrators
	PDESWEPlaneMoriZwanzigTimeSteppers pdeSWEPlaneMoriZwanzigTimeSteppers;

	// Handler to all benchmarks
	PDESWEPlaneMoriZwanzigBenchmarksCombined planeBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWEPlaneMoriZwanzig *shackPDESWEPlane;
	ShackPDESWEPlaneMoriZwanzigTimeDiscretization *shackTimeDisc;
	ShackPDESWEPlaneMoriZwanzigBenchmarks *shackPDESWEPlaneBenchmarks;
	ShackPDESWEPlaneDiagnostics *shackPDESWEPlaneDiagnostics;

	class BenchmarkErrors
	{
	public:
		//Max difference to initial conditions
		double t0_diff_error_max_abs_h_pert_SP;
		double t0_diff_error_max_abs_u_SP;
		double t0_diff_error_max_abs_v_SP;
		double t0_diff_error_max_abs_h_pert_SQ;
		double t0_diff_error_max_abs_u_SQ;
		double t0_diff_error_max_abs_v_SQ;
		double t0_diff_error_max_abs_h_pert_FQ;
		double t0_diff_error_max_abs_u_FQ;
		double t0_diff_error_max_abs_v_FQ;

		// Error measures L2 norm
		double analytical_error_rms_h_SP;
		double analytical_error_rms_u_SP;
		double analytical_error_rms_v_SP;
		double analytical_error_rms_h_SQ;
		double analytical_error_rms_u_SQ;
		double analytical_error_rms_v_SQ;
		double analytical_error_rms_h_FQ;
		double analytical_error_rms_u_FQ;
		double analytical_error_rms_v_FQ;

		// Error measures max norm
		double analytical_error_maxabs_h_SP;
		double analytical_error_maxabs_u_SP;
		double analytical_error_maxabs_v_SP;
		double analytical_error_maxabs_h_SQ;
		double analytical_error_maxabs_u_SQ;
		double analytical_error_maxabs_v_SQ;
		double analytical_error_maxabs_h_FQ;
		double analytical_error_maxabs_u_FQ;
		double analytical_error_maxabs_v_FQ;

		void setup()
		{
			t0_diff_error_max_abs_h_pert_SP = -1;
			t0_diff_error_max_abs_u_SP = -1;
			t0_diff_error_max_abs_v_SP = -1;
			t0_diff_error_max_abs_h_pert_SQ = -1;
			t0_diff_error_max_abs_u_SQ = -1;
			t0_diff_error_max_abs_v_SQ = -1;
			t0_diff_error_max_abs_h_pert_FQ = -1;
			t0_diff_error_max_abs_u_FQ = -1;
			t0_diff_error_max_abs_v_FQ = -1;

			analytical_error_rms_h_SP = -1;
			analytical_error_rms_u_SP = -1;
			analytical_error_rms_v_SP = -1;
			analytical_error_rms_h_SQ = -1;
			analytical_error_rms_u_SQ = -1;
			analytical_error_rms_v_SQ = -1;
			analytical_error_rms_h_FQ = -1;
			analytical_error_rms_u_FQ = -1;
			analytical_error_rms_v_FQ = -1;

			analytical_error_maxabs_h_SP = -1;
			analytical_error_maxabs_u_SP = -1;
			analytical_error_maxabs_v_SP = -1;
			analytical_error_maxabs_h_SQ = -1;
			analytical_error_maxabs_u_SQ = -1;
			analytical_error_maxabs_v_SQ = -1;
			analytical_error_maxabs_h_FQ = -1;
			analytical_error_maxabs_u_FQ = -1;
			analytical_error_maxabs_v_FQ = -1;
		}
	};

	BenchmarkErrors benchmarkErrors;

#if 0
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
#endif

	/// Diagnostic measures at initial stage, Initialize with 0
	double diagnostics_energy_start = 0;
	double diagnostics_mass_start = 0;
	double diagnostics_potential_enstrophy_start = 0;


	bool compute_error_difference_to_initial_condition = false;
	bool compute_error_to_analytical_solution = false;
	bool compute_normal_modes = false;


public:
	ProgramPDESWEPlaneMoriZwanzig(
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
		shackPDESWEPlane = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzig>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzigTimeDiscretization>();
		shackPDESWEPlaneBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneMoriZwanzigBenchmarks>();
		shackPDESWEPlaneDiagnostics = shackProgArgDict.getAutoRegistration<ShackPDESWEPlaneDiagnostics>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		planeBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);

		pdeSWEPlaneMoriZwanzigTimeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWEPlaneMoriZwanzigTimeSteppers);

#if 0
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes = new NormalModesData;
			normalmodes->shackRegistration(shackProgArgDict);
		}
#endif

		return true;
	}

	void clear_1_shackRegistration()
	{
#if 0
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			delete normalmodes;
			normalmodes = nullptr;
		}
#endif

		shackPlaneDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		planeBenchmarksCombined.clear();
		pdeSWEPlaneMoriZwanzigTimeSteppers.clear();
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
		dataAndOps_SP.setup(shackPlaneDataOps);
		dataAndOps_SQ.setup(shackPlaneDataOps);
		dataAndOps_FQ.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataAndOps_SP);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataAndOps_SQ);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataAndOps_FQ);

		/*
		 * After we setup the plane, we can setup the time steppers and their buffers
		 */
		pdeSWEPlaneMoriZwanzigTimeSteppers.setup(&dataAndOps_SP.ops, &shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(pdeSWEPlaneMoriZwanzigTimeSteppers);

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		planeBenchmarksCombined.setupInitialConditions(
				dataAndOps_SP.prog_h_pert,
				dataAndOps_SP.prog_u,
				dataAndOps_SP.prog_v,
				&dataAndOps_SP.ops,
				&dataAndOps_SP.planeDataConfig
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);


		// copy and project initial solutions
		dataAndOps_SQ.prog_h_pert = dataAndOps_SP.prog_h_pert;
		dataAndOps_SQ.prog_u = dataAndOps_SP.prog_u;
		dataAndOps_SQ.prog_v = dataAndOps_SP.prog_v;
		dataAndOps_FQ.prog_h_pert = dataAndOps_SP.prog_h_pert;
		dataAndOps_FQ.prog_u = dataAndOps_SP.prog_u;
		dataAndOps_FQ.prog_v = dataAndOps_SP.prog_v;

		projection.project_SP(dataAndOps_SP.prog_h_pert, dataAndOps_SP.prog_u, dataAndOps_SP.prog_v);
		projection.project_SQ(dataAndOps_SQ.prog_h_pert, dataAndOps_SQ.prog_u, dataAndOps_SQ.prog_v);
		projection.project_FQ(dataAndOps_FQ.prog_h_pert, dataAndOps_FQ.prog_u, dataAndOps_FQ.prog_v);


		dataAndOps_SP.t0_prog_h_pert = dataAndOps_SP.prog_h_pert;
		dataAndOps_SQ.t0_prog_h_pert = dataAndOps_SQ.prog_h_pert;
		dataAndOps_FQ.t0_prog_h_pert = dataAndOps_FQ.prog_h_pert;

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

#if 0
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes->setup(&dataAndOps.planeDataConfig);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*normalmodes);
		}
#endif

		if (shackPDESWEPlane->compute_errors)
		{
			//Compute difference to initial condition (makes more sense in steady state cases, but useful in others too)
			compute_error_difference_to_initial_condition = true;

			//Compute difference to analytical solution (makes more sense in linear cases, but might be useful in others too)
			compute_error_to_analytical_solution = pdeSWEPlaneMoriZwanzigTimeSteppers.linear_only;
		}
		else
		{
			compute_error_difference_to_initial_condition = false;
			compute_error_to_analytical_solution = false;
		}



#if 0
		/*
		 * Load initial conditions from file if required
		 */
		if (shackIOData->initial_condition_data_filenames.size() > 0)
			dataAndOps.prog_h_pert.file_physical_loadData(shackIOData->initial_condition_data_filenames[0].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 1)
			dataAndOps.prog_u.file_physical_loadData(shackIOData->initial_condition_data_filenames[1].c_str(), shackIOData->initial_condition_input_data_binary);

		if (shackIOData->initial_condition_data_filenames.size() > 2)
			dataAndOps.prog_v.file_physical_loadData(shackIOData->initial_condition_data_filenames[2].c_str(), shackIOData->initial_condition_input_data_binary);
#endif

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
#if 0
		if (shackPDESWEPlane->normal_mode_analysis_generation)
		{
			normalmodes->clear();
			delete normalmodes;
			normalmodes = nullptr;
		}
#endif
		
#if SWEET_GUI
		vis_plane_data.clear();
#endif

		pdeSWEPlaneMoriZwanzigTimeSteppers.clear();

		dataAndOps_SP.clear();
		dataAndOps_SQ.clear();
		dataAndOps_FQ.clear();
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

#if 0
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

#endif
	}


	//Update diagnostic variables related to normal modes
	void dump_normal_modes()
	{
#if 0
		if (!compute_normal_modes )
			return;

#if SWEET_USE_PLANE_SPECTRAL_SPACE
		PDESWEPlaneBench_NormalModes n;
		n.dump_all_normal_modes(normalmodes->geo, normalmodes->igwest, normalmodes->igeast);
#endif

		return;
#endif
	}

	//Calculate the model diagnostics
	void update_diagnostics()
	{

#if 0
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
#endif
	}



	void normal_mode_analysis()
	{
#if 0
		normalmodes->pdeSWEPlaneNormalModes.normal_mode_analysis(
								dataAndOps.prog_h_pert,
								dataAndOps.prog_u,
								dataAndOps.prog_v,
								3,
								&shackProgArgDict,
								this,
								&ProgramPDESWEPlaneMoriZwanzig::runTimestep
						);

		std::cout << "\n Done normal mode analysis in separate class" << std::endl;
#endif
	}


	/**
	 * Execute a single simulation time step
	 */
	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();

		pdeSWEPlaneMoriZwanzigTimeSteppers.timestepper->runTimestep(
				dataAndOps_SP.prog_h_pert, dataAndOps_SP.prog_u, dataAndOps_SP.prog_v,
				dataAndOps_SQ.prog_h_pert, dataAndOps_SQ.prog_u, dataAndOps_SQ.prog_v,
				dataAndOps_FQ.prog_h_pert, dataAndOps_FQ.prog_u, dataAndOps_FQ.prog_v,
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
		sweet::PlaneData_Physical t_h(dataAndOps_SP.planeDataConfig);
		sweet::PlaneData_Physical t_u(dataAndOps_SP.planeDataConfig);
		sweet::PlaneData_Physical t_v(dataAndOps_SP.planeDataConfig);

		SimDataAndOps* dataAndOps;

		std::vector<std::string> sol_cases = {"SP", "SQ", "FQ"};
		for (std::vector<std::string>::iterator it = sol_cases.begin(); it != sol_cases.end(); it++)
		{
			std::string sol_case = *it;

			if (sol_case == "SP")
				dataAndOps = &dataAndOps_SP;
			else if (sol_case == "SQ")
				dataAndOps = &dataAndOps_SQ;
			else if (sol_case == "FQ")
				dataAndOps = &dataAndOps_FQ;

			if (shackPlaneDataOps->space_grid_use_c_staggering) // Remap in case of C-grid
			{
				t_h = dataAndOps->prog_h_pert.toPhys();
				dataAndOps->gridMapping.mapCtoA_u(dataAndOps->prog_u.toPhys(), t_u);
				dataAndOps->gridMapping.mapCtoA_v(dataAndOps->prog_v.toPhys(), t_v);
			}
			else
			{
				t_h = dataAndOps->prog_h_pert.toPhys();
				t_u = dataAndOps->prog_u.toPhys();
				t_v = dataAndOps->prog_v.toPhys();
			}

			update_normal_modes();

			// Dump  data in csv, if output filename is not empty
			if (shackIOData->output_file_name.size() > 0)
			{
				output_filenames = "";

				output_filenames = write_file(t_h, ("prog_h_pert_" + sol_case).c_str());
				output_filenames += ";" + write_file(t_u, ("prog_u_" + sol_case).c_str());
				output_filenames += ";" + write_file(t_v, ("prog_v_" + sol_case).c_str());

				output_filenames += ";" + write_file(dataAndOps->ops.ke(t_u,t_v), ("diag_ke_" + sol_case).c_str());

#if SWEET_USE_PLANE_SPECTRAL_SPACE
				output_filenames += ";" + write_file_spec(dataAndOps->ops.ke(t_u,t_v), ("diag_ke_" + sol_case + "_spec").c_str());

				output_filenames += ";" + write_file_spec(t_h, ("prog_h_pert_" + sol_case + "_spec").c_str());
				output_filenames += ";" + write_file_spec(t_u, ("prog_u_" + sol_case + "_spec").c_str());
				output_filenames += ";" + write_file_spec(t_v, ("prog_v_" + sol_case + "_spec").c_str());

				output_filenames += ";" + write_file_spec(dataAndOps->ops.ke(t_u,t_v).toPhys(), ("diag_ke_" + sol_case + "_spec").c_str());
#endif

				output_filenames += ";" + write_file(dataAndOps->ops.vort(t_u, t_v), ("diag_vort_" + sol_case).c_str());
				output_filenames += ";" + write_file(dataAndOps->ops.div(t_u, t_v), ("diag_div_" + sol_case).c_str());

////////////
////////#if SWEET_USE_PLANE_SPECTRAL_SPACE
////////				if (compute_normal_modes){
////////					output_filenames += ";" + write_file_spec(normalmodes->geo, "nm_geo");
////////					output_filenames += ";" + write_file_spec(normalmodes->igwest, "nm_igwest");
////////					output_filenames += ";" + write_file_spec(normalmodes->igeast, "nm_igeast");
////////				}
////////#endif
			}

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

#if 0
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

#if 0
			if (compute_error_difference_to_initial_condition)
			{
				// Difference to initial condition
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tDIFF_MAXABS_H0\tDIFF_MAXABS_U0\tDIFF_MAXABS_V0";

				rows << "\t" << benchmarkErrors.t0_diff_error_max_abs_h_pert << "\t" << benchmarkErrors.t0_diff_error_max_abs_u << "\t" << benchmarkErrors.t0_diff_error_max_abs_v;
			}
#endif

#if 0
			if (compute_error_to_analytical_solution)
			{
				if (shackTimestepControl->current_timestep_nr == 0)
					header << "\tREF_DIFF_MAX_H\tREF_DIFF_MAX_U\tREF_DIFF_MAX_V";

				rows << "\t" << benchmarkErrors.analytical_error_maxabs_h << "\t" << benchmarkErrors.analytical_error_maxabs_u << "\t" << benchmarkErrors.analytical_error_maxabs_v;
			}
#endif

#if 0
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

#if 0

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

#endif

	}


public:
	void printErrors()
	{

#if 0
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

#endif

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

	bool instability_detected()
	{
		return !(	
					dataAndOps_SP.prog_h_pert.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_SP.prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_SP.prog_v.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_SQ.prog_h_pert.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_SQ.prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_SQ.prog_v.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_FQ.prog_h_pert.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_FQ.prog_u.toPhys().physical_reduce_boolean_all_finite() &&
					dataAndOps_FQ.prog_v.toPhys().physical_reduce_boolean_all_finite()
				);
	}

};

#endif
