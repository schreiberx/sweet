/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_SWE_SPHERE_HPP_
#define SRC_PROGRAMS_PDE_SWE_SPHERE_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/sphere/Sphere.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "ShackPDESWESphere.hpp"

// Benchmarks
#include "PDESWESphereBenchmarksCombined.hpp"

// Time steppers
#include "PDESWESphereTimeSteppers.hpp"

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
	#include <sweet/core/plane/Plane.hpp>
	#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#endif

#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/sphere/SphereData_DebugContainer.hpp>




class ProgramPDESWESphere
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

		sweet::SphereData_Config sphereDataConfig;
		sweet::SphereOperators ops;

		sweet::SphereData_Spectral prog_phi_pert;
		sweet::SphereData_Spectral prog_div;
		sweet::SphereData_Spectral prog_vrt;


		sweet::SphereData_Spectral t0_prog_phi_pert;
		sweet::SphereData_Spectral t0_prog_div;
		sweet::SphereData_Spectral t0_prog_vrt;

#if SWEET_GUI
	sweet::PlaneData_Config planeDataConfig;

	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical vis_plane_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 1;

	// Which primitive to use for rendering
	int vis_data_id = 0;

#endif


		bool setup(
				sweet::ShackSphereDataOps *i_shackSphereDataOps
		)
		{
			/*
			 * Setup Sphere Data Config & Operators
			 */
			sphereDataConfig.setupAuto(i_shackSphereDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(sphereDataConfig);

			ops.setup(&sphereDataConfig, i_shackSphereDataOps);
			ERROR_CHECK_WITH_RETURN_BOOLEAN(ops);

			prog_phi_pert.setup(sphereDataConfig);
			prog_div.setup(sphereDataConfig);
			prog_vrt.setup(sphereDataConfig);

#if SWEET_GUI
			sweet::ShackPlaneDataOps shackPlaneDataOps;
			shackPlaneDataOps.space_res_physical[0] = i_shackSphereDataOps->space_res_physical[0];
			shackPlaneDataOps.space_res_physical[1] = i_shackSphereDataOps->space_res_physical[1];
			shackPlaneDataOps.reuse_spectral_transformation_plans = i_shackSphereDataOps->reuse_spectral_transformation_plans;

			planeDataConfig.setupAuto(shackPlaneDataOps);
#endif
			return true;
		}

		void clear()
		{
			prog_phi_pert.clear();
			prog_div.clear();
			prog_vrt.clear();

			t0_prog_phi_pert.clear();
			t0_prog_div.clear();
			t0_prog_vrt.clear();

			ops.clear();
			sphereDataConfig.clear();
		}
	};

	// Simulation data
	Data dataConfigOps;

	// time integrators
	PDESWESphereTimeSteppers timeSteppers;

	// Handler to all benchmarks
	PDESWESphereBenchmarksCombined benchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWESphereTimeDiscretization *shackTimeDisc;
	ShackPDESWESphereBenchmarks *shackBenchmarks;
	ShackPDESWESphere *shackPDESWESphere;



public:
	ProgramPDESWESphere(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackBenchmarks(nullptr),
		shackPDESWESphere(nullptr)
	{
		ERROR_CHECK_WITH_RETURN(shackProgArgDict);
	}


	bool setup_1_shackRegistration()
	{
		/*
		 * Setup argument parsing
		 */
		shackProgArgDict.setup();

		/*
		 * SHACK: Register classes which we require
		 */
		shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWESphereBenchmarks>();
		shackPDESWESphere = shackProgArgDict.getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */

		/*
		 * Setup benchmarks
		 */
		benchmarksCombined.setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(benchmarksCombined);

		benchmarksCombined.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Setup time steppers
		 */
		timeSteppers.setup_1_registerAllTimesteppers();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(timeSteppers);

		timeSteppers.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(timeSteppers);

		/*
		 * Process HELP arguments
		 */
		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Close shack registration & getting shacks
		 */
		shackProgArgDict.closeRegistration();
		shackProgArgDict.closeGet();

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphereDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		benchmarksCombined.clear();
		timeSteppers.clear();
		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}

	bool setup_3_data()
	{
		/*
		 * BENCHMARK: Detect particular benchmark to use
		 */
		benchmarksCombined.setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(benchmarksCombined);

		/*
		 * Setup benchmark itself
		 */
		benchmarksCombined.setup_4_benchmarkSetup_1_withoutOps();

		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphereDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(dataConfigOps);

		/*
		 * Setup benchmark itself
		 */
		benchmarksCombined.setup_5_benchmarkSetup_2_withOps(&dataConfigOps.ops);

		/*
		 * Now we're ready to setup the time steppers
		 */
		timeSteppers.setup_3_timestepper(
				shackTimeDisc->timestepping_method,
				&shackProgArgDict,
				&dataConfigOps.ops
			);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(timeSteppers);


		/*
		 * Load initial state of benchmark
		 */
		benchmarksCombined.benchmark->getInitialState(
				dataConfigOps.prog_phi_pert,
				dataConfigOps.prog_vrt,
				dataConfigOps.prog_div
			);

		ERROR_CHECK_WITH_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Backup data at t=0
		 */
		dataConfigOps.t0_prog_phi_pert = dataConfigOps.prog_phi_pert;
		dataConfigOps.t0_prog_vrt = dataConfigOps.prog_vrt;
		dataConfigOps.t0_prog_div = dataConfigOps.prog_div;

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

		return true;
	}
	void clear_3_data()
	{
#if SWEET_GUI
		dataConfigOps.vis_plane_data.clear();
#endif

		timeSteppers.clear();

		dataConfigOps.clear();
	}

	bool setup()
	{
		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_data())
			return false;

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		std::cout << "SETUP FINISHED" << std::endl;
		return true;
	}
	void clear()
	{
		clear_3_data();
		clear_2_processArguments();
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

	void printSimulationErrors()
	{
		sweet::SphereData_Spectral diff = dataConfigOps.t0_prog_phi_pert-dataConfigOps.prog_phi_pert;

		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << diff.toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << "RMS error: " << diff.toPhys().physical_reduce_rms() << std::endl;
	}

	~ProgramPDESWESphere()
	{
		clear();
	}


	bool runTimestep()
	{
		sweet::SphereData_DebugContainer::clear();

		shackTimestepControl->timestepHelperStart();


		timeSteppers.timestepper->runTimestep(
				dataConfigOps.prog_phi_pert, dataConfigOps.prog_vrt, dataConfigOps.prog_div,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);

#if 0
		if (shackBenchmarks->getVelocities)
		{
			/*
			 * From advection test case
			 *
			 * Update velocities just for sake of the correction visualization
			 */
			shackBenchmarks->getVelocities(
					dataConfigOps.vel_u,
					dataConfigOps.vel_v,
					shackTimestepControl->current_simulation_time + shackTimestepControl->current_timestep_size,
					shackBenchmarks->getVelocitiesUserData
				);
		}
#endif

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 2)
			std::cout << shackTimestepControl->current_timestep_nr << ": " << shackTimestepControl->current_simulation_time/(60*60*24.0) << std::endl;

		return true;
	}



	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}

};







////////////////////////////////////////////////
#ifndef SWEET_GUI
	#define SWEET_GUI 1
#endif


#include <stdexcept>

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
	#include <sweet/core/plane/PlaneData_Config.hpp>
	#include <sweet/core/plane/PlaneData_Physical.hpp>
	#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
	#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>
#endif

#include "PDESWESphereBenchmarksCombined.hpp"

#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereData_Physical.hpp>
#include <sweet/core/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/sphere/SphereOperatorsComplex.hpp>
#include <sweet/core/sphere/SphereData_SpectralComplex.hpp>

#include <sweet/core/Stopwatch.hpp>
#include <sweet/core/SWEETError.hpp>

#include "PDESWESphereTimeSteppers.hpp"
#include "PDESWESphereNormalModeAnalysis.hpp"

#include <sweet/core/StopwatchBox.hpp>
#include <sweet/core/sphere/SphereData_DebugContainer.hpp>

#if SWEET_PARAREAL
	#include <sweet/parareal/Parareal.hpp>
#endif

#if SWEET_XBRAID
	#include <sweet/xbraid/XBraid_sweet_lib.hpp>
#endif

sweet::ShackDictionary shackDict;

// Plane data config
sweet::SphereData_Config sphereDataConfigInstance;
sweet::SphereData_Config sphereDataConfigInstance_nodealiasing;
sweet::SphereData_Config *sphereDataConfig = &sphereDataConfigInstance;
sweet::SphereData_Config *sphereDataConfig_nodealiasing = &sphereDataConfigInstance_nodealiasing;


#if SWEET_GUI
	sweet::PlaneData_Config planeDataConfigInstance;
	sweet::PlaneData_Config *planeDataConfig = &planeDataConfigInstance;
#endif



/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class ProgramPDESWESphereXXX
{
public:
	sweet::SphereOperators op;
	sweet::SphereOperators op_nodealiasing;

	PDESWESphereTimeSteppers timeSteppers;

#if SWEET_PARAREAL
	// Implementation of different time steppers
	PDESWESphereTimeSteppers timeSteppersCoarse;
#endif


	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;

	sweet::SphereData_Spectral prog_phi_pert;
	sweet::SphereData_Spectral prog_vrt;
	sweet::SphereData_Spectral prog_div;

	Stopwatch stopwatch;

#if SWEET_GUI
	sweet::PlaneData_Physical viz_plane_data;
#endif

	int render_primitive_id = 1;

	sweet::SphereHelpers_Diagnostics sphereDiagnostics;

#if SWEET_MPI
	int mpi_rank;
#endif

	// was the output of the time step already done for this simulation state?
	double timestep_last_output_simtime;

	PDESWESphereBenchmarksCombined sphereBenchmarks;

public:
	ProgramPDESWESphere]()	:
		op(sphereDataConfig, &(shackDict.sim)),
		op_nodealiasing(sphereDataConfig_nodealiasing, &(shackDict.sim)),
		prog_phi_pert(sphereDataConfig),
		prog_vrt(sphereDataConfig),
		prog_div(sphereDataConfig),

#if SWEET_GUI
		viz_plane_data(planeDataConfig),
#endif
		sphereDiagnostics(
				sphereDataConfig,
				shackDict,
				shackDict.misc.verbosity
		)

	{
#if SWEET_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif


		reset();

	}



	void update_diagnostics()
	{
		// assure, that the diagnostics are only updated for new time steps
		if (last_timestep_nr_update_diagnostics == shackDict.timecontrol.current_timestep_nr)
			return;

		sphereDiagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
				op,
				prog_phi_pert,
				prog_vrt,
				prog_div,
				shackDict
		);
	}



	void reset()
	{
		StopwatchBox::getInstance().main_setup.start();

		shackDict.reset();
		shackDict.iodata.output_time_scale = 1.0/(60.0*60.0);

		// Diagnostics measures
		last_timestep_nr_update_diagnostics = -1;

		shackDict.iodata.output_next_sim_seconds = 0;

		if (shackDict.timecontrol.current_timestep_size <= 0)
			SWEETError("Only fixed time step size supported");

		if (shackDict.benchmark.setup_dealiased)
		{
			//std::cout << "A" << std::endl;
			//std::cout << sphereDataConfig->getConfigInformationString() << std::endl;
			//exit(1);
			// use dealiased physical space for setup
			sphereBenchmarks.setup(shackDict, op);
			sphereBenchmarks.timestepper->getInitialState(prog_phi_pert, prog_vrt, prog_div);
		}
		else
		{
			//std::cout << "B" << std::endl;
			//std::cout << sphereDataConfig_nodealiasing->getConfigInformationString() << std::endl;
			//exit(1);
			// this is not the default since noone uses it
			// use reduced physical space for setup to avoid spurious modes
			sweet::SphereData_Spectral prog_phi_pert_nodealiasing(sphereDataConfig_nodealiasing);
			sweet::SphereData_Spectral prog_vrt_nodealiasing(sphereDataConfig_nodealiasing);
			sweet::SphereData_Spectral prog_div_nodealiasing(sphereDataConfig_nodealiasing);

			sphereBenchmarks.setup(shackDict, op_nodealiasing);
			sphereBenchmarks.timestepper->getInitialState(prog_phi_pert_nodealiasing, prog_vrt_nodealiasing, prog_div_nodealiasing);

			prog_phi_pert.load_nodealiasing(prog_phi_pert_nodealiasing);
			prog_vrt.load_nodealiasing(prog_vrt_nodealiasing);
			prog_div.load_nodealiasing(prog_div_nodealiasing);
		}

		/*
		 * SETUP time steppers
		 */
		timeSteppers.setup(shackDict.disc.timestepping_method,
				op, shackDict);

		std::cout << "[MULE] timestepper_string_id: " << timeSteppers.timestepper->string_id() << std::endl;

		update_diagnostics();

		shackDict.diag.backup_reference();

		StopwatchBox::getInstance().main_setup.stop();

		// start at one second in the past to ensure output at t=0
		timestep_last_output_simtime = shackDict.timecontrol.current_simulation_time-1.0;

		/*
		 * Output configuration here to ensure that updated variables are included in this output
		 */
#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			shackDict.outputConfig();
		}


		/*
		 * Output data for the first time step as well if output of datafiels is requested
		 */
		if (shackDict.iodata.output_each_sim_seconds >= 0)
			timestep_do_output();

		stopwatch.start();
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv_spec_evol(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name		///< name of output variable
	)
	{
		char buffer[1024];
		std::string phase = "_arg";
		std::string ampl = "_amp"; 

		const char* filename_template_ampl = "output_spec_ampl_%s.txt"; //.c_str();
		const char* filename_template_arg = "output_spec_arg_%s.txt"; //.c_str();
		int reduce_mode_factor = 4;

		sprintf(buffer, filename_template_arg, i_name);
		i_sphereData.spectrum_phase_file_write_line(buffer, 
			i_name, shackDict.timecontrol.current_simulation_time*shackDict.iodata.output_time_scale,
			20, 10e-20, reduce_mode_factor);

		sprintf(buffer, filename_template_ampl, i_name);
		i_sphereData.spectrum_abs_file_write_line(buffer, 
			i_name, shackDict.timecontrol.current_simulation_time*shackDict.iodata.output_time_scale,
			20, 10e-20, reduce_mode_factor);

		return buffer;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name,		///< name of output variable
			bool i_phi_shifted = false
	)
	{
		char buffer[1024];

		// create copy
		sweet::SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = shackDict.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackDict.timecontrol.current_simulation_time*shackDict.iodata.output_time_scale);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name
	)
	{
		char buffer[1024];

		sweet::SphereData_Spectral sphereData(i_sphereData);
		const char* filename_template = shackDict.iodata.output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackDict.timecontrol.current_simulation_time*shackDict.iodata.output_time_scale);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return;
#endif

		if (shackDict.iodata.output_file_name.length() == 0)
			return;

		std::cout << "Writing output files at simulation time: " << shackDict.timecontrol.current_simulation_time << " secs" << std::endl;

		if (shackDict.iodata.output_file_mode == "csv")
		{
			std::string output_filename;

			sweet::SphereData_Spectral h = prog_phi_pert*(1.0/shackDict.sim.gravitation);
			h += shackDict.sim.h0;

			output_filename = write_file_csv(h, "prog_h");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;

			output_filename = write_file_csv(prog_phi_pert, "prog_phi_pert");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << prog_phi_pert.toPhys().physical_reduce_min() << ", max: " << prog_phi_pert.toPhys().physical_reduce_max() << ")" << std::endl;

			sweet::SphereData_Physical phi_phys = h.toPhys() * shackDict.sim.gravitation;
			sweet::SphereData_Spectral phi(sphereDataConfig);
			phi.loadSphereDataPhysical(phi_phys);
			output_filename = write_file_csv(phi, "prog_phi");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << phi_phys.physical_reduce_min() << ", max: " << phi_phys.physical_reduce_max() << ")" << std::endl;

			sweet::SphereData_Physical u(sphereDataConfig);
			sweet::SphereData_Physical v(sphereDataConfig);

			op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);

			output_filename = write_file_csv(u, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(v, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(prog_vrt, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(prog_div, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			sweet::SphereData_Spectral potvrt = (prog_phi_pert/shackDict.sim.gravitation)*prog_vrt;

			output_filename = write_file_csv(potvrt, "prog_potvrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;
		}
		else if (shackDict.iodata.output_file_mode == "bin")
		{
			std::string output_filename;

			{
				output_filename = write_file_bin(prog_phi_pert, "prog_phi_pert");
				output_reference_filenames = output_filename;
				sweet::SphereData_Physical prog_phys = prog_phi_pert.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(prog_vrt, "prog_vrt");
				output_reference_filenames += ";"+output_filename;
				sweet::SphereData_Physical prog_phys = prog_vrt.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(prog_div, "prog_div");
				output_reference_filenames += ";"+output_filename;
				sweet::SphereData_Physical prog_phys = prog_div.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}
		}
		else if (shackDict.iodata.output_file_mode == "csv_spec_evol"){

			std::string output_filename;

			{ 
				/*
				* Spectral kinetic energy and potential enstrophy calculation and output
				*
				* Details in Jakob-Chien, Ruediger, James J. Hack, and David L. Williamson. 
				* "Spectral transform solutions to the shallow water test set." Journal of Computational Physics 119, no. 1 (1995): 164-187.
				*/
				// Kinetic energy is given in spectral space as
				// KE per mode = a^2/((n(n+1)))*(vrt*conj(vrt))+a^2/((n(n+1)))*(div*conj(div))
				// r = a/(sqrt(n(n+1))) (root_laplace)
				// KE per mode = (r*vrt*conj(r*vrt))+(r*div*conj(r*div))
				sweet::SphereData_Spectral rlap_vrt = op.inv_root_laplace(prog_vrt); 
				sweet::SphereData_Spectral rlap_div = op.inv_root_laplace(prog_div); 
				sweet::SphereData_Spectral kin_en = rlap_vrt + rlap_div ;

				output_filename = write_file_csv_spec_evol(kin_en, "kin_en"); 
				std::cout << " + " << output_filename << " (Total Kin Energy : " << 0.25*kin_en.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// For Barotropic vort eq: See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// Kinetic energy is given in spectral space as
				// Vortical energy per mode is (0.5 n*(n+1) / a^2) *psi*conj(psi) in spectral space
				//SphereData_Spectral psi = op.inv_laplace(prog_vrt); // 
				// multiply psi by sqrt( n * (n+1))/a (apply root laplacian)
				//SphereData_Spectral psi_root = op.root_laplace(psi);
				//output_filename = write_file_csv_spec_evol(psi_root*std::sqrt(0.5), "spec_energy"); 
				//std::cout << " + " << output_filename << " (Kinetic energy : " << (0.5)*psi_root.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// enstrophy per mode is 0.5 vrt*conj(vrt) in spectral space
				// Total enstrophy is the sum of these (counting twice modes with m>0 and once when m=0)
				output_filename = write_file_csv_spec_evol(prog_vrt, "enstrophy"); 
				std::cout << " + " << output_filename << " (Total Enstrophy : " << prog_vrt.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

			}
		}
		else
		{
			SWEETError("Unknown output file mode '"+shackDict.iodata.output_file_mode+"'");
		}
	}



	void timestep_do_output()
	{
		if (shackDict.misc.compute_errors)
		{
			/*
			 * Check for stationary solutions
			 */
			if (
					shackDict.benchmark.benchmark_name != "williamson2"					&&
					shackDict.benchmark.benchmark_name != "williamson2_linear"			&&
					shackDict.benchmark.benchmark_name != "galewsky_nobump"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_linear"	&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_1"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_2"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_4"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_8"			&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_16"		&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_32"		&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_64"		&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_128"		&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_256"		&&
					shackDict.benchmark.benchmark_name != "geostrophic_balance_512"
			)
			{
				std::cout << "Benchmark name: " << shackDict.benchmark.benchmark_name << std::endl;
				SWEETError("Analytical solution not available for this benchmark");
			}

			sweet::SphereData_Spectral anal_solution_phi_pert(sphereDataConfig);
			sweet::SphereData_Spectral anal_solution_vrt(sphereDataConfig);
			sweet::SphereData_Spectral anal_solution_div(sphereDataConfig);

			sphereBenchmarks.setup(shackDict, op);
			sphereBenchmarks.timestepper->getInitialState(anal_solution_phi_pert, anal_solution_vrt, anal_solution_div);

			/*
			 * Compute difference
			 */
			sweet::SphereData_Spectral diff_phi = prog_phi_pert - anal_solution_phi_pert;
			sweet::SphereData_Spectral diff_vrt = prog_vrt - anal_solution_vrt;
			sweet::SphereData_Spectral diff_div = prog_div - anal_solution_div;

#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				double error_phi = diff_phi.toPhys().physical_reduce_max_abs();
				double error_vrt = diff_vrt.toPhys().physical_reduce_max_abs();
				double error_div = diff_div.toPhys().physical_reduce_max_abs();

				
				std::ios init(NULL);
				init.copyfmt(std::cout);
				std::cout << "[MULE] errors." << std::setw(8) << std::setfill('0') << shackDict.timecontrol.current_timestep_nr << ": ";
				std::cout.copyfmt(init);

				std::cout << "simtime=" << shackDict.timecontrol.current_simulation_time;
				std::cout << "\terror_linf_phi=" << error_phi;
				std::cout << "\terror_linf_vrt=" << error_vrt;
				std::cout << "\terror_linf_div=" << error_div;
				std::cout << std::endl;
			}
		}

		write_file_output();

		update_diagnostics();

		if (shackDict.misc.verbosity > 1)
		{

#if SWEET_MPI
			if (mpi_rank == 0)
#endif
			{
				update_diagnostics();

				// Print header
				if (shackDict.timecontrol.current_timestep_nr == 0)
				{
					std::cout << "T\tTOTAL_MASS\tPOT_ENERGY\tKIN_ENERGY\tTOT_ENERGY\tPOT_ENSTROPHY\tREL_TOTAL_MASS\tREL_POT_ENERGY\tREL_KIN_ENERGY\tREL_TOT_ENERGY\tREL_POT_ENSTROPHY";
					std::cout << std::endl;
				}

				// Print simulation time, energy and pot enstrophy
				std::cout << shackDict.timecontrol.current_simulation_time << "\t";
				std::cout << shackDict.diag.total_mass << "\t";
				std::cout << shackDict.diag.potential_energy << "\t";
				std::cout << shackDict.diag.kinetic_energy << "\t";
				std::cout << shackDict.diag.total_energy << "\t";
				std::cout << shackDict.diag.total_potential_enstrophy << "\t";

				std::cout << (shackDict.diag.total_mass-shackDict.diag.ref_total_mass)/shackDict.diag.total_mass << "\t";
				std::cout << (shackDict.diag.potential_energy-shackDict.diag.ref_potential_energy)/shackDict.diag.potential_energy << "\t";
				std::cout << (shackDict.diag.kinetic_energy-shackDict.diag.ref_kinetic_energy)/shackDict.diag.kinetic_energy << "\t";
				std::cout << (shackDict.diag.total_energy-shackDict.diag.total_energy)/shackDict.diag.total_energy << "\t";
				std::cout << (shackDict.diag.total_potential_enstrophy-shackDict.diag.total_potential_enstrophy)/shackDict.diag.total_potential_enstrophy << std::endl;

				static double start_tot_energy = -1;
				if (start_tot_energy == -1)
					start_tot_energy = shackDict.diag.total_energy;
			}
		}


		if (shackDict.misc.verbosity > 0)
		{
#if SWEET_MPI
			if (mpi_rank == 0)
#endif
				std::cout << "prog_phi min/max:\t" << prog_phi_pert.toPhys().physical_reduce_min() << ", " << prog_phi_pert.toPhys().physical_reduce_max() << std::endl;
		}

		if (shackDict.iodata.output_each_sim_seconds > 0)
			while (shackDict.iodata.output_next_sim_seconds <= shackDict.timecontrol.current_simulation_time)
				shackDict.iodata.output_next_sim_seconds += shackDict.iodata.output_each_sim_seconds;
	}




public:
	bool timestep_check_output()
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return false;
#endif
		if (shackDict.misc.gui_enabled)
			update_diagnostics();

		if (shackDict.misc.verbosity > 0)
			std::cout << "." << std::flush;

		if (shackDict.iodata.output_each_sim_seconds < 0)
			return false;

		if (shackDict.timecontrol.current_simulation_time == timestep_last_output_simtime)
			return false;

		timestep_last_output_simtime = shackDict.timecontrol.current_simulation_time;

		if (shackDict.timecontrol.current_simulation_time < shackDict.timecontrol.max_simulation_time - shackDict.iodata.output_each_sim_seconds*1e-10)
		{
			if (shackDict.iodata.output_next_sim_seconds > shackDict.timecontrol.current_simulation_time)
				return false;
		}

		if (shackDict.misc.verbosity > 0)
			std::cout << std::endl;

		timestep_do_output();

		return true;
	}



public:
	bool should_quit()
	{
		if (shackDict.timecontrol.max_timesteps_nr != -1 && shackDict.timecontrol.max_timesteps_nr <= shackDict.timecontrol.current_timestep_nr)
			return true;

		if (shackDict.timecontrol.max_wallclock_time >= 0)
		{
			double t = stopwatch.getTimeSinceStart();
			if (shackDict.timecontrol.max_wallclock_time <= t)
			{
				std::cout << "[MULE] max_wallclock_time: exceeded (" << t << ")" << std::endl;
				return true;
			}
		}

		double diff = std::abs(shackDict.timecontrol.max_simulation_time - shackDict.timecontrol.current_simulation_time);

		if (	shackDict.timecontrol.max_simulation_time != -1 &&
				(
						shackDict.timecontrol.max_simulation_time <= shackDict.timecontrol.current_simulation_time	||
						diff/shackDict.timecontrol.max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
				)
			)
			return true;

		return false;
	}



	bool detect_instability()
	{
		if (prog_phi_pert.spectral_is_first_nan_or_inf())
		{
			std::cout << "Infinity value detected" << std::endl;
			std::cerr << "Infinity value detected" << std::endl;
			return true;
		}

		return false;
	}



	void runTimestep()
	{
#if SWEET_GUI
		if (shackDict.misc.gui_enabled && shackDict.misc.normal_mode_analysis_generation == 0)
			timestep_check_output();
#endif

		if (shackDict.timecontrol.current_simulation_time + shackDict.timecontrol.current_timestep_size > shackDict.timecontrol.max_simulation_time)
			shackDict.timecontrol.current_timestep_size = shackDict.timecontrol.max_simulation_time - shackDict.timecontrol.current_simulation_time;

		timeSteppers.timestepper->runTimestep(
				prog_phi_pert, prog_vrt, prog_div,
				shackDict.timecontrol.current_timestep_size,
				shackDict.timecontrol.current_simulation_time
			);



		// Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral

		if (shackDict.sim.viscosity != 0 && shackDict.misc.use_nonlinear_only_visc == 0)
		{
			///prog_vrt = op.implicit_diffusion(prog_vrt, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.sphere_radius);
			///prog_div = op.implicit_diffusion(prog_div, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.sphere_radius);
			///prog_phi_pert = op.implicit_diffusion(prog_phi_pert, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.sphere_radius);
			prog_vrt = op.implicit_hyperdiffusion(prog_vrt, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.viscosity_order, shackDict.sim.sphere_radius);
			prog_div = op.implicit_hyperdiffusion(prog_div, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.viscosity_order, shackDict.sim.sphere_radius);
			prog_phi_pert = op.implicit_hyperdiffusion(prog_phi_pert, shackDict.timecontrol.current_timestep_size*shackDict.sim.viscosity, shackDict.sim.viscosity_order, shackDict.sim.sphere_radius);
		}


		// advance time step and provide information to parameters
		shackDict.timecontrol.current_simulation_time += shackDict.timecontrol.current_timestep_size;
		shackDict.timecontrol.current_timestep_nr++;

#if SWEET_GUI
		timestep_check_output();
#endif
	}



	void normalmode_analysis()
	{
		NormalModeAnalysisSphere::normal_mode_analysis(
				prog_phi_pert,
				prog_vrt,
				prog_div,
				shackDict,
				this,
				&ProgramPDESWESphere]::run_timestep
			);
	}



#if SWEET_GUI

	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		if (shackDict.timecontrol.run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				runTimestep();
	}


	int max_viz_types = 9;


	void vis_get_vis_data_array(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		// request rendering of sphere
		*o_render_primitive_id = render_primitive_id;
		*o_bogus_data = sphereDataConfig;

		if (shackDict.misc.vis_id < 0)
		{
			int n = -shackDict.misc.vis_id-1;
			if (n <  (int)SphereData_DebugContainer().size())
			{
				SphereData_DebugContainer::DataContainer &d = sweet::SphereData_DebugContainer().container_data()[n];
				if (d.is_spectral)
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(d.data_spectral, planeDataConfig);
				else
					viz_plane_data = Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(d.data_physical, planeDataConfig);

				*o_dataArray = &viz_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
		}

		int id = shackDict.misc.vis_id % max_viz_types;

		switch (id)
		{
			default:

			case 0:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_phi_pert), planeDataConfig);
				break;

			case 1:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_vrt), planeDataConfig);
				break;

			case 2:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(SphereData_Spectral(prog_div), planeDataConfig);
				break;

			case 3:
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(shackDict.sim.h0 + sweet::SphereData_Spectral(prog_phi_pert)/shackDict.sim.gravitation, planeDataConfig);
				break;

			case 4:
			{
				sweet::SphereData_Physical u(prog_vrt.sphereDataConfig);
				sweet::SphereData_Physical v(prog_vrt.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(u, planeDataConfig);
				break;
			}

			case 5:
			{
				sweet::SphereData_Physical u(prog_vrt.sphereDataConfig);
				sweet::SphereData_Physical v(prog_vrt.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				op.vrtdiv_to_uv(prog_vrt, prog_div, u, v);
				viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(v, planeDataConfig);
				break;
			}

			case 6:
			case 7:
			case 8:
			{
				sweet::SphereData_Spectral anal_solution_phi_pert(sphereDataConfig);
				sweet::SphereData_Spectral anal_solution_vrt(sphereDataConfig);
				sweet::SphereData_Spectral anal_solution_div(sphereDataConfig);

				sphereBenchmarks.setup(shackDict, op);
				sphereBenchmarks.timestepper->getInitialState(anal_solution_phi_pert, anal_solution_vrt, anal_solution_div);

				switch (id)
				{
				case 6:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_phi_pert - anal_solution_phi_pert, planeDataConfig);
					break;

				case 7:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_vrt - anal_solution_vrt, planeDataConfig);
					break;

				case 8:
					viz_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(prog_div - anal_solution_div, planeDataConfig);
					break;
				}
			}
		}


		double viz_min = viz_plane_data.physical_reduce_min();
		double viz_max = viz_plane_data.physical_reduce_max();

		viz_max = std::max(std::abs(viz_max), std::abs(viz_min));
		viz_min = -viz_max;

		*o_viz_min = viz_min;
		*o_viz_max = viz_max;


		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const char* vis_getStatusString()
	{
		std::string description = "";


		bool found = false;
		if (shackDict.misc.vis_id < 0)
		{
			int n = -shackDict.misc.vis_id-1;

			if (n <  (int)SphereData_DebugContainer().size())
			{
				description = std::string("DEBUG_")+SphereData_DebugContainer().container_data()[n].description;
				found = true;
			}
		}

		int id = shackDict.misc.vis_id % max_viz_types;

		if (!found)
		{
			switch (id)
			{
			default:
			case 0:
				description = "phi_pert";
				break;

			case 1:
				description = "vrt";
				break;

			case 2:
				description = "div";
				break;

			case 3:
				description = "h";
				break;

			case 4:
				description = "u";
				break;

			case 5:
				description = "v";
				break;

			case 6:
				description = "phi diff t0";
				break;

			case 7:
				description = "vrt diff t0";
				break;

			case 8:
				description = "div diff t0";
				break;
			}
		}


		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
#if SWEET_MPI
				"Rank %i,"
				","
#endif
				"Visualization %i: %s,"
				"MaxVal: %.12e,"
				"MinVal: %.12e,"
				","
				"Time: %f secs,"
				"Time: %f hours,"
				"Time: %f days,"
				"timestep nr.: %i,"
				"timestep size: %f,"
				","
				"TMass: %.12e,"
				"TEnergy: %.12e,"
				"PotEnstrophy: %.12e,"
				","
				"Colorscale: lowest [Blue... green ... red] highest"
				,
#if SWEET_MPI
				mpi_rank,
#endif
				shackDict.misc.vis_id,
				description.c_str(),
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min(),

				shackDict.timecontrol.current_simulation_time,
				shackDict.timecontrol.current_simulation_time/(60.0*60.0),
				shackDict.timecontrol.current_simulation_time/(60.0*60.0*24.0),
				shackDict.timecontrol.current_timestep_nr,
				shackDict.timecontrol.current_timestep_size,

				shackDict.diag.total_mass,
				shackDict.diag.total_energy,
				shackDict.diag.total_potential_enstrophy

		);

		return title_string;
	}



	void vis_pause()
	{
		shackDict.timecontrol.run_simulation_timesteps = !shackDict.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			shackDict.misc.vis_id++;
			break;

		case 'V':
			shackDict.misc.vis_id--;
			break;

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
			break;

		case 'c':
			write_file_output();
			break;
		}
	}
#endif


};


int main_real(int i_argc, char *i_argv[])
{

	// Time counter
	StopwatchBox::getInstance().main.start();

#if __MIC__
	std::cout << "Compiled for MIC" << std::endl;
#endif

#if SWEET_MPI

	#if SWEET_THREADING_SPACE
		int provided;
		MPI_Init_thread(&i_argc, &i_argv, MPI_THREAD_MULTIPLE, &provided);

		if (provided != MPI_THREAD_MULTIPLE)
		{
				std::cerr << "MPI_THREAD_MULTIPLE not available! Try to get an MPI version with multi-threading support or compile without OMP/TBB support. Good bye..." << std::endl;
				exit(-1);
		}
	#else
		MPI_Init(&i_argc, &i_argv);
	#endif

#endif

	// input parameter names (specific ones for this program)
	const char *user_defined_prog_params[] = {
			nullptr
	};

#if SWEET_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

	// Help menu
	if (!shackDict.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params))
	{
#if SWEET_PARAREAL
		shackDict.parareal.printProgramArguments();
#endif
		return -1;
	}

	if (shackDict.misc.verbosity > 3)
		std::cout << " + setup SH sphere transformations..." << std::endl;

	sphereDataConfigInstance.setupAuto(
			shackDict.disc.space_res_physical,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackDict.misc.verbosity,
			shackDict.parallelization.num_threads_space
		);

	int res_physical_nodealias[2] = {
			2*shackDict.disc.space_res_spectral[0],
			shackDict.disc.space_res_spectral[1]
		};

	if (shackDict.misc.verbosity > 3)
		std::cout << " + setup SH sphere transformations (nodealiasing)..." << std::endl;

	sphereDataConfigInstance_nodealiasing.setupAuto(
			res_physical_nodealias,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackDict.misc.verbosity,
			shackDict.parallelization.num_threads_space
		);


#if SWEET_GUI
	if (shackDict.misc.verbosity > 3)
		std::cout << " + setup FFT plane transformations..." << std::endl;

	planeDataConfigInstance.setupAutoSpectralSpaceFromPhysical(shackDict.disc.space_res_physical, shackDict.misc.reuse_spectral_transformation_plans);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);

	if (shackDict.misc.verbosity > 3)
		std::cout << " + setup finished" << std::endl;

#if SWEET_MPI
	std::cout << "Hello from MPI rank: " << mpi_rank << std::endl;

	// only start simulation and time stepping for first rank
	if (mpi_rank > 0)
	{
		/*
		 * Deactivate all output for ranks larger than the current one
		 */
		shackDict.misc.verbosity = 0;
	#if !SWEET_XBRAID
		shackDict.iodata.output_each_sim_seconds = -1;
	#endif
	}
#endif

	{
#if SWEET_MPI
		if (mpi_rank == 0)
#endif
		{
			std::cout << "SPH config string: " << sphereDataConfigInstance.getConfigInformationString() << std::endl;
		}

#if SWEET_PARAREAL
		if (shackDict.parareal.enabled)
		{

			shackDict.iodata.output_time_scale = 1.0/(60.0*60.0);

			//SphereOperators op(sphereDataConfig, shackDict.sim.plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs);
			sweet::SphereOperators op(sphereDataConfig, &(shackDict.sim));
			sweet::SphereOperators op_nodealiasing(sphereDataConfig_nodealiasing, &(shackDict.sim));

			// Set planeDataConfig and planeOperators for each level
			std::vector<SphereData_Config*> sphereDataConfigs;
			std::vector<SphereOperators*> ops;
			std::vector<SphereOperators*> ops_nodealiasing;

			// fine
			sphereDataConfigs.push_back(sphereDataConfig);
			ops.push_back(&op);
			ops_nodealiasing.push_back(&op_nodealiasing);

			// coarse
			if (shackDict.parareal.spatial_coarsening)
			{
				///for (int j = 0; j < 2; j++)
				///	assert(shackDict.disc.space_res_physical[j] == -1);
				int N_physical[2] = {-1, -1};
				int N_spectral[2];

				double frac;
				if ( shackDict.parareal.coarse_timestep_size > 0)
					frac = shackDict.timecontrol.current_timestep_size / shackDict.parareal.coarse_timestep_size;
				else
					frac = shackDict.timecontrol.current_timestep_size / (shackDict.timecontrol.max_simulation_time / shackDict.parareal.coarse_slices );
				for (int j = 0; j < 2; j++)
					N_spectral[j] = std::max(4, int(shackDict.disc.space_res_spectral[j] * frac));

				sphereDataConfigs.push_back(new sweet::SphereData_Config);
				sphereDataConfigs.back()->setupAuto(
						N_physical,
						N_spectral,
						shackDict.misc.reuse_spectral_transformation_plans,
						shackDict.misc.verbosity,
						shackDict.parallelization.num_threads_space
					);

				ops.push_back(new sweet::SphereOperators(sphereDataConfigs.back(), &(shackDict.sim)));
				// @TODO: nodealiasing case
				ops_nodealiasing.push_back(ops.back());
			}
			else
			{
				sphereDataConfigs.push_back(sphereDataConfig);
				ops.push_back(&op);
				ops_nodealiasing.push_back(&op_nodealiasing);
			}


			PDESWESphereTimeSteppers* timeSteppersFine = new PDESWESphereTimeSteppers;
			PDESWESphereTimeSteppers* timeSteppersCoarse = new PDESWESphereTimeSteppers;

			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller<PDESWESphereTimeSteppers, 3> parareal_Controller(	&shackDict,
												sphereDataConfigs,
												ops,
												ops_nodealiasing,
												timeSteppersFine,
												timeSteppersCoarse);

			// setup controller. This initializes several simulation instances
			parareal_Controller.setup();

			// execute the simulation
			parareal_Controller.run();

			delete timeSteppersFine;
			delete timeSteppersCoarse;

			if (shackDict.parareal.spatial_coarsening)
			{
				delete sphereDataConfigs[1];
				delete ops[1];
				sphereDataConfigs[1] = nullptr;
				ops[1] = nullptr;
			}

		}
		else
#endif

#if SWEET_XBRAID

		if (shackDict.xbraid.xbraid_enabled)
		{

			shackDict.iodata.output_time_scale = 1.0/(60.0*60.0);

			sweet::SphereOperators op(sphereDataConfig, &(shackDict.sim));

			// Set planeDataConfig and planeOperators for each level
			std::vector<SphereData_Config*> sphereDataConfigs;
			std::vector<SphereOperators*> ops;
			for (int i = 0; i < shackDict.xbraid.xbraid_max_levels; i++)
			{
				if (shackDict.xbraid.xbraid_spatial_coarsening)
				{
					int N_physical[2] = {-1, -1};
					int N_spectral[2];
					for (int j = 0; j < 2; j++)
					{
						// proportional to time step
						if (shackDict.xbraid.xbraid_spatial_coarsening == 1)
							N_spectral[j] = std::max(4,int(shackDict.disc.space_res_spectral[j] / std::pow(shackDict.xbraid.xbraid_cfactor, i)));
						else if (shackDict.xbraid.xbraid_spatial_coarsening > 1)
						{
							if (i == 0)
								N_spectral[j] = std::max(4, shackDict.disc.space_res_spectral[j]);
							else
								N_spectral[j] = std::max(4, shackDict.xbraid.xbraid_spatial_coarsening);
						}
						else
							SWEETError("Invalid parameter xbraid_spatial_coarsening");
					}

					sphereDataConfigs.push_back(new sweet::SphereData_Config);
					sphereDataConfigs.back()->setupAuto(
							N_physical,
							N_spectral,
							shackDict.misc.reuse_spectral_transformation_plans,
							shackDict.misc.verbosity,
							shackDict.parallelization.num_threads_space
						);

					ops.push_back(new sweet::SphereOperators(sphereDataConfigs.back(), &(shackDict.sim)));

					std::cout << "Spectral resolution at level " << i << " : " << N_spectral[0] << " " << N_spectral[1] << std::endl;
				}
				else
				{
					sphereDataConfigs.push_back(sphereDataConfig);
					ops.push_back(&op);
				}
			}






			MPI_Comm comm = MPI_COMM_WORLD;
			MPI_Comm comm_x, comm_t;

			//////braid_Core core;
			///sweet_App* app = (sweet_App *) malloc(sizeof(sweet_App))
			int nt = (int) (shackDict.timecontrol.max_simulation_time / shackDict.timecontrol.current_timestep_size);
                        if (nt * shackDict.timecontrol.current_timestep_size < shackDict.timecontrol.max_simulation_time - 1e-10)
				nt++;
			///sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackDict.timecontrol.max_simulation_time, nt, &shackDict, sphereDataConfig, &op);
			sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackDict.timecontrol.max_simulation_time, nt, &shackDict, sphereDataConfigs, ops);


			if ( shackDict.xbraid.xbraid_run_wrapper_tests)
			{
				app.setup();

				BraidUtil braid_util;
				int test = braid_util.TestAll(&app, comm, stdout, 0., shackDict.timecontrol.current_timestep_size, shackDict.timecontrol.current_timestep_size * 2);
				////int test = braid_util.TestBuf(app, comm, stdout, 0.);
				if (test == 0)
					SWEETError("Tests failed!");
				else
					std::cout << "Tests successful!" << std::endl;

			}
			else
			{
				BraidCore core(MPI_COMM_WORLD, &app);
				app.setup(core);
				// Run Simulation
				core.Drive();
			}

			if (shackDict.xbraid.xbraid_spatial_coarsening)
				for (int i = 0; i < shackDict.xbraid.xbraid_max_levels; i++)
				{
					delete sphereDataConfigs[i];
					delete ops[i];
					sphereDataConfigs[i] = nullptr;
					ops[i] = nullptr;
				}

		}
		else
#endif


#if SWEET_GUI // The VisSweet directly calls simulationSWE->reset() and output stuff
		if (shackDict.misc.gui_enabled)
		{
			ProgramPDESWESphere] *simulationSWE = new ProgramPDESWESphere];
			VisSweet<ProgramPDESWESphere]> visSweet(simulationSWE);
			delete simulationSWE;
		}
		else
#endif
		{
			ProgramPDESWESphere] *simulationSWE = new ProgramPDESWESphere];

			if (shackDict.misc.normal_mode_analysis_generation > 0)
			{
				simulationSWE->normalmode_analysis();
			}
			else
			{
				// Do first output before starting timer
				simulationSWE->timestep_check_output();
#if SWEET_MPI
				// Start counting time
				if (mpi_rank == 0)
				{
					std::cout << "********************************************************************************" << std::endl;
					std::cout << "Parallel performance information: MPI barrier & timer starts here" << std::endl;
					std::cout << "********************************************************************************" << std::endl;
				}
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				StopwatchBox::getInstance().main_timestepping.start();

				// Main time loop
				while (true)
				{
					// Stop simulation if requested
					if (simulationSWE->should_quit())
						break;

					// Test for some output to be done
					simulationSWE->timestep_check_output();

					// Main call for timestep run
					simulationSWE->runTimestep();

					// Instability
					if (shackDict.misc.instability_checks)
					{
#if SWEET_MPI
						if (mpi_rank == 0)
#endif
						{
							if (simulationSWE->detect_instability())
							{
								std::cout << "INSTABILITY DETECTED" << std::endl;
								std::cerr << "INSTABILITY DETECTED" << std::endl;
								// IMPORANT: EXIT IN CASE OF INSTABILITIES
								exit(1);
								break;
							}
						}
					}
				}

				// Stop counting time
				StopwatchBox::getInstance().main_timestepping.stop();

#if SWEET_MPI
				MPI_Barrier(MPI_COMM_WORLD);
#endif

				if (shackDict.misc.verbosity > 0)
					std::cout << std::endl;
#if SWEET_MPI
				// Start counting time
				if (mpi_rank == 0)
				{
					std::cout << "********************************************************************************" << std::endl;
					std::cout << "Parallel performance information: timer stopped here" << std::endl;
					std::cout << "********************************************************************************" << std::endl;
				}
#endif

				// Do some output after the time loop
				simulationSWE->timestep_check_output();
			}

#if SWEET_MPI
			// Start counting time
			if (mpi_rank == 0)
#endif
			{
				if (shackDict.iodata.output_file_name.size() > 0)
					std::cout << "[MULE] reference_filenames: " << simulationSWE->output_reference_filenames << std::endl;
			}

			std::cout << "[MULE] simulation_successfully_finished: 1" << std::endl;

			delete simulationSWE;
		}

		StopwatchBox::getInstance().main.stop();
	}


#if SWEET_MPI
	if (mpi_rank == 0)
#endif
	{
		// End of run output results
		std::cout << std::endl;
		StopwatchBox::getInstance().output();

		std::cout << "***************************************************" << std::endl;
		std::cout << "* Other timing information (direct)" << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "[MULE] shackDict.timecontrol.current_timestep_nr: " << shackDict.timecontrol.current_timestep_nr << std::endl;
		std::cout << "[MULE] shackDict.timecontrol.current_timestep_size: " << shackDict.timecontrol.current_timestep_size << std::endl;
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "* Other timing information (derived)" << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "[MULE] simulation_benchmark_timings.time_per_time_step (secs/ts): " << StopwatchBox::getInstance().main_timestepping()/(double)shackDict.timecontrol.current_timestep_nr << std::endl;
	}

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}

#endif
