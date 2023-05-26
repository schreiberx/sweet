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
#include <sweet/core/shacksShared/ShackParallelization.hpp>
#include "benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "ShackPDESWESphere.hpp"


#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
	#include <sweet/core/plane/Plane.hpp>
	#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#endif

#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/StopwatchBox.hpp>

// Benchmarks
#include "PDESWESphere_BenchmarksCombined.hpp"

// Time steppers
#include "time/PDESWESphere_TimeSteppers.hpp"
// Time steppers
#include "timeTree/PDESWESphere_TimeSteppersNewTS.hpp"

#include "timeTree/PDESWESphere_DataContainer.hpp"

// Diagnostics
#include "PDESWESphere_Diagnostics.hpp"

// Normal mode analysis
#include "PDESWESphere_NormalModeAnalysis.hpp"

// File writer
#include "PDESWESphere_FileOutput.hpp"


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
	class DataConfigOps
	{
	public:
		sweet::ErrorBase error;

		sweet::SphereData_Config sphereDataConfig;
		sweet::SphereOperators ops;

		PDESWESphere_DataContainer prog;
		PDESWESphere_DataContainer progTmp;

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
				sweet::ShackSphereDataOps *i_shackSphereDataOps,
				bool i_setup_spectral_transforms = true		// for reset()
		)
		{
			/*
			 * Setup Sphere Data Config & Operators
			 */
			if (i_setup_spectral_transforms)
			{
				sphereDataConfig.setupAuto(i_shackSphereDataOps);
				ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereDataConfig);
			}

			ops.setup(&sphereDataConfig, i_shackSphereDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog.setup(&sphereDataConfig);
			progTmp.setup(&sphereDataConfig);

#if SWEET_GUI
			sweet::ShackPlaneDataOps shackPlaneDataOps;
			shackPlaneDataOps.space_res_physical[0] = i_shackSphereDataOps->space_res_physical[0];
			shackPlaneDataOps.space_res_physical[1] = i_shackSphereDataOps->space_res_physical[1];
			shackPlaneDataOps.reuse_spectral_transformation_plans = i_shackSphereDataOps->reuse_spectral_transformation_plans;

			planeDataConfig.setupAuto(shackPlaneDataOps);
#endif
			return true;
		}

		void clear(bool i_clear_spectral_transforms = true)
		{
			prog.clear();

			t0_prog_phi_pert.clear();
			t0_prog_div.clear();
			t0_prog_vrt.clear();

			ops.clear();

			if (i_clear_spectral_transforms)
				sphereDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// time integrators
	PDESWESphere_TimeSteppers timeSteppers;
	PDESWESphere_TimeSteppersNewTS timeSteppersNewTS;
	bool useNewTimeSteppers;

	// Handler to all benchmarks
	PDESWESphere_BenchmarksCombined sphereBenchmarks;

	PDESWESphere_FileOutput fileOutput;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	sweet::ShackParallelization *shackParallelization;
	ShackPDESWESphereTimeDiscretization *shackTimeDisc;
	ShackPDESWESphereBenchmarks *shackBenchmarks;
	ShackPDESWESphere *shackPDESWESphere;
	
	int timestep_nr_last_output_simtime = -1;

	PDESWESphere_Diagnostics diagnostics;



public:
	ProgramPDESWESphere(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackParallelization(nullptr),
		shackTimeDisc(nullptr),
		shackBenchmarks(nullptr),
		shackPDESWESphere(nullptr)
	{
		ERROR_CHECK_COND_RETURN(shackProgArgDict);
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
		shackParallelization = shackProgArgDict.getAutoRegistration<sweet::ShackParallelization>();
		shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDESWESphereBenchmarks>();
		shackPDESWESphere = shackProgArgDict.getAutoRegistration<ShackPDESWESphere>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */

		/*
		 * Setup benchmarks
		 */
		sphereBenchmarks.setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereBenchmarks);

		sphereBenchmarks.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereBenchmarks);

		shackProgArgDict.processProgramArguments(shackTimeDisc);


		if (shackTimeDisc->timestepping_method.find('(') != std::string::npos)
			useNewTimeSteppers = true;
		else
			useNewTimeSteppers = false;

		if (useNewTimeSteppers)
		{
			/*
			 * Setup NEW time steppers
			 */
			timeSteppersNewTS.setup_1_registerAllTimesteppers();
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);
		}
		else
		{
			/*
			 * Setup legacy time steppers
			 */
			timeSteppers.setup_1_registerAllTimesteppers();
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

			timeSteppers.setup_2_shackRegistration(&shackProgArgDict);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);
		}

		/*
		 * Process HELP arguments
		 */
		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		if (!useNewTimeSteppers)
		{
			/*
			 * Close shack registration & getting shacks
			 */
			shackProgArgDict.closeRegistration();

			shackProgArgDict.closeGet();
		}

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackSphereDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackParallelization = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		sphereBenchmarks.clear();
		timeSteppers.clear();
		shackProgArgDict.clear();
	}


	bool setup_2_processArguments()
	{
		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * Do some validation of program arguments
		 */
		shackTimestepControl->validateTimestepSize();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(*shackTimestepControl);

		return true;
	}

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}


	bool setup_3_dataAndOps(bool i_setup_spectral_transforms = true)
	{
		/*
		 * BENCHMARK: Detect particular benchmark to use
		 */
		sphereBenchmarks.setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereBenchmarks);

		/*
		 * Setup benchmark itself
		 */
		sphereBenchmarks.setup_4_benchmarkSetup_1_withoutOps();

		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphereDataOps, i_setup_spectral_transforms);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		/*
		 * Setup benchmark itself
		 */
		sphereBenchmarks.setup_5_benchmarkSetup_2_withOps(&dataConfigOps.ops);

		/*
		 * Now we're ready to setup the time steppers
		 */
		if (useNewTimeSteppers)
		{
			timeSteppersNewTS.setup_2_timestepper(
				shackTimeDisc->timestepping_method,
				&shackProgArgDict,
				&dataConfigOps.ops,
				dataConfigOps.prog
			);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);

			timeSteppersNewTS.timeIntegrator->setTimeStepSize(shackTimestepControl->current_timestep_size);

			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppersNewTS);
		}
		else
		{
			timeSteppers.setup_3_timestepper(
					shackTimeDisc->timestepping_method,
					&shackProgArgDict,
					&dataConfigOps.ops
				);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);
		}
		
		if (useNewTimeSteppers)
		{
			shackProgArgDict.closeRegistration();
			shackProgArgDict.closeGet();
		}


		/*
		 * Load initial state of benchmark
		 */
		sphereBenchmarks.benchmark->getInitialState(
				dataConfigOps.prog.phi_pert,
				dataConfigOps.prog.vrt,
				dataConfigOps.prog.div
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereBenchmarks);


		/*
		 * Backup data at t=0
		 */
		dataConfigOps.t0_prog_phi_pert = dataConfigOps.prog.phi_pert;
		dataConfigOps.t0_prog_vrt = dataConfigOps.prog.vrt;
		dataConfigOps.t0_prog_div = dataConfigOps.prog.div;

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
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		if (shackPDESWESphere->compute_diagnostics)
			diagnostics.setup(&(dataConfigOps.ops), shackPDESWESphere, 0);

		fileOutput.setup(shackIOData, shackTimestepControl, shackPDESWESphere);

		return true;
	}
	void clear_3_data(bool i_clear_spectral_transforms = true)
	{
		fileOutput.clear();

#if SWEET_GUI
		dataConfigOps.vis_plane_data.clear();
#endif

		if (useNewTimeSteppers)
		{
			timeSteppersNewTS.clear();
		}
		else
		{
			timeSteppers.clear();
		}

		dataConfigOps.clear(i_clear_spectral_transforms);
	}

	bool setup(bool i_setup_spectral_transforms = true)
	{
		StopwatchBox::getInstance().main_setup.start();

		if (!setup_1_shackRegistration())
			return false;

		if (!setup_2_processArguments())
			return false;

		if (!setup_3_dataAndOps(i_setup_spectral_transforms))
			return false;

		if (shackParallelization->isMPIRoot)
		{
			std::cout << "Printing shack information:" << std::endl;

			shackProgArgDict.printShackData();

			std::cout << "SETUP FINISHED" << std::endl;
		}

		StopwatchBox::getInstance().main_setup.stop();

		/*
		 * Output data for the first time step as well if output of datafiels is requested
		 */
		if (shackIOData->output_each_sim_seconds >= 0)
			_timestepDoOutput();
		return true;
	}
	void clear(bool i_clear_spectral_transforms = true)
	{
		clear_3_data(i_clear_spectral_transforms);
		clear_2_processArguments();
		clear_1_shackRegistration();
	}

	bool reset()
	{
		// keep pausing simulation
		bool run_simulation_timesteps = shackTimestepControl->run_simulation_timesteps;

		clear(false);

		if (!setup(false))
		{
			error.print();
			return false;
		}

		shackTimestepControl->run_simulation_timesteps = run_simulation_timesteps;

		return !error.exists();
	}

	void printSimulationErrors()
	{
		sweet::SphereData_Spectral diff = dataConfigOps.t0_prog_phi_pert-dataConfigOps.prog.phi_pert;

		double lmax_error = diff.toPhys().physical_reduce_max_abs();
		double rms_error = diff.toPhys().physical_reduce_rms();

		if (shackParallelization->isMPIRoot)
		{
			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << lmax_error << std::endl;
			std::cout << "RMS error: " << rms_error << std::endl;
		}
	}

	~ProgramPDESWESphere()
	{
		clear();
	}


	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();


		if (useNewTimeSteppers)
		{
			timeSteppersNewTS.timeIntegrator->eval_integration(
					dataConfigOps.prog,
					dataConfigOps.progTmp,
					shackTimestepControl->current_simulation_time
				);
			dataConfigOps.prog.swap(dataConfigOps.progTmp);
		}
		else
		{
			timeSteppers.timestepper->runTimestep(
					dataConfigOps.prog.phi_pert, dataConfigOps.prog.vrt, dataConfigOps.prog.div,
					shackTimestepControl->current_timestep_size,
					shackTimestepControl->current_simulation_time
				);
		}

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 2)
			if (shackParallelization->isMPIRoot)
				std::cout << shackTimestepControl->current_timestep_nr << ": " << shackTimestepControl->current_simulation_time/(60*60*24.0) << std::endl;

		if (shackPDESWESphere->compute_diagnostics)
			update_diagnostics();

		return true;
	}



	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}



	void update_diagnostics()
	{
		if (diagnostics.last_update_timestep_nr == shackTimestepControl->current_timestep_nr)
			return;

		assert(shackPDESWESphere->compute_diagnostics);

		diagnostics.update_phi_vrt_div_2_mass_energy_enstrophy(
				&dataConfigOps.ops,
				dataConfigOps.prog.phi_pert,
				dataConfigOps.prog.vrt,
				dataConfigOps.prog.div,

				shackSphereDataOps->sphere_radius,
				shackPDESWESphere->gravitation
		);
	}


	void _timestepDoOutput()
	{
		if (shackPDESWESphere->compute_diagnostics)
		{
			update_diagnostics();

			if (shackParallelization->isMPIRoot)
			{
				// Print header
				if (shackTimestepControl->current_timestep_nr == 0)
					diagnostics.printTabularHeader();

				diagnostics.printTabularRow(shackTimestepControl->current_simulation_time);
			}
		}

		if (shackPDESWESphere->compute_errors)
		{
			/*
			 * Check for stationary solutions
			 */
			if (
				shackBenchmarks->benchmark_name != "williamson2"				&&
				shackBenchmarks->benchmark_name != "williamson2_linear"			&&
				shackBenchmarks->benchmark_name != "galewsky_nobump"			&&
				shackBenchmarks->benchmark_name != "geostrophic_balance"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_linear"	&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_1"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_2"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_4"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_8"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_16"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_32"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_64"		&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_128"	&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_256"	&&
				shackBenchmarks->benchmark_name != "geostrophic_balance_512"
			)
			{

				if (shackParallelization->isMPIRoot)
				{
					std::cout << "Benchmark name: " << shackBenchmarks->benchmark_name << std::endl;
				}
				SWEETError("Analytical solution not available for this benchmark");
			}

			sweet::SphereData_Spectral anal_solution_phi_pert(dataConfigOps.sphereDataConfig);
			sweet::SphereData_Spectral anal_solution_vrt(dataConfigOps.sphereDataConfig);
			sweet::SphereData_Spectral anal_solution_div(dataConfigOps.sphereDataConfig);

			sphereBenchmarks.benchmark->getInitialState(anal_solution_phi_pert, anal_solution_vrt, anal_solution_div);

			/*
			 * Compute difference
			 */
			sweet::SphereData_Spectral diff_phi = dataConfigOps.prog.phi_pert - anal_solution_phi_pert;
			sweet::SphereData_Spectral diff_vrt = dataConfigOps.prog.vrt - anal_solution_vrt;
			sweet::SphereData_Spectral diff_div = dataConfigOps.prog.div - anal_solution_div;

			double error_phi = diff_phi.toPhys().physical_reduce_max_abs();
			double error_vrt = diff_vrt.toPhys().physical_reduce_max_abs();
			double error_div = diff_div.toPhys().physical_reduce_max_abs();

			if (shackParallelization->isMPIRoot)
			{
				int nTimeSteps = shackTimestepControl->current_timestep_nr;
				std::cout << "[MULE] errors." << std::setw(8) << std::setfill('0') << nTimeSteps << ": ";

				std::cout << "simtime=" << shackTimestepControl->current_simulation_time;
				std::cout << "\terror_linf_phi=" << error_phi;
				std::cout << "\terror_linf_vrt=" << error_vrt;
				std::cout << "\terror_linf_div=" << error_div;
				std::cout << std::endl;
			}
		}

		if (shackParallelization->isMPIRoot)
		{
			fileOutput.write_file_output(
					dataConfigOps.ops,
					dataConfigOps.prog.phi_pert,
					dataConfigOps.prog.div,
					dataConfigOps.prog.vrt
			);

			if (shackIOData->verbosity > 0)
			{
				double progPhiMin = dataConfigOps.prog.phi_pert.toPhys().physical_reduce_min();
				double progPhiMax = dataConfigOps.prog.phi_pert.toPhys().physical_reduce_max();

				if (shackParallelization->isMPIRoot)
				{
					std::cout << "prog_phi min/max:\t" << progPhiMin << ", " << progPhiMax << std::endl;
				}
			}
		}

		if (shackIOData->output_each_sim_seconds > 0)
			while (shackIOData->output_next_sim_seconds <= shackTimestepControl->current_simulation_time)
				shackIOData->output_next_sim_seconds += shackIOData->output_each_sim_seconds;
	}



public:
	bool timestepHandleOutput()
	{
		if (shackIOData->output_each_sim_seconds < 0)
			return false;

		if (shackTimestepControl->current_simulation_time == timestep_nr_last_output_simtime)
			return false;

		timestep_nr_last_output_simtime = shackTimestepControl->current_simulation_time;

		if (shackTimestepControl->current_simulation_time < shackTimestepControl->max_simulation_time - shackIOData->output_each_sim_seconds*1e-10)
		{
			if (shackIOData->output_next_sim_seconds > shackTimestepControl->current_simulation_time)
				return false;
		}

		if (shackParallelization->isMPIRoot)
			if (shackIOData->verbosity > 0)
				std::cout << std::endl;

		_timestepDoOutput();

		return true;
	}



	bool detect_instability()
	{
		if (dataConfigOps.prog.phi_pert.spectral_is_first_nan_or_inf())
		{
			if (shackParallelization->isMPIRoot)
			{
				std::cout << "Infinity value detected" << std::endl;
				std::cerr << "Infinity value detected" << std::endl;
			}
			return true;
		}

		return false;
	}

	void normalmode_analysis()
	{
		if (useNewTimeSteppers)
		{
			SWEETError("Not supported, yet");
		}
		else
		{
		PDESWESphere_NormalModeAnalysis::normal_mode_analysis(
				dataConfigOps.prog.phi_pert,
				dataConfigOps.prog.vrt,
				dataConfigOps.prog.div,

				shackIOData,
				shackTimestepControl,

				shackPDESWESphere->normal_mode_analysis_generation,

				this,
				&ProgramPDESWESphere::runTimestep
			);
		}
	}

	void output_timings()
	{
		if (shackParallelization->isMPIRoot)
		{
			std::cout << std::endl;
			StopwatchBox::getInstance().output();

			std::cout << "***************************************************" << std::endl;
			std::cout << "* Other timing information (direct)" << std::endl;
			std::cout << "***************************************************" << std::endl;
			std::cout << "[MULE] shackTimestepControl->current_timestep_nr: " << shackTimestepControl->current_timestep_nr << std::endl;
			std::cout << "[MULE] shackTimestepControl->current_timestep_size: " << shackTimestepControl->current_timestep_size << std::endl;
			std::cout << std::endl;
			std::cout << "***************************************************" << std::endl;
			std::cout << "* Other timing information (derived)" << std::endl;
			std::cout << "***************************************************" << std::endl;
			std::cout << "[MULE] simulation_benchmark_timings.time_per_time_step (secs/ts): " << StopwatchBox::getInstance().main_timestepping()/(double)shackTimestepControl->current_timestep_nr << std::endl;
		}
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


	int max_vis_types = 9;


	void vis_getDataArray(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_vis_render_type_of_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		// request rendering of sphere
		*o_vis_render_type_of_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
		*o_bogus_data = &dataConfigOps.sphereDataConfig;

		int id = dataConfigOps.vis_data_id % max_vis_types;

		switch (id)
		{
			default:

			case 0:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(sweet::SphereData_Spectral(dataConfigOps.prog.phi_pert), dataConfigOps.planeDataConfig);
				break;

			case 1:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(sweet::SphereData_Spectral(dataConfigOps.prog.vrt), dataConfigOps.planeDataConfig);
				break;

			case 2:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(sweet::SphereData_Spectral(dataConfigOps.prog.div), dataConfigOps.planeDataConfig);
				break;

			case 3:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(shackPDESWESphere->h0 + sweet::SphereData_Spectral(dataConfigOps.prog.phi_pert)/shackPDESWESphere->gravitation, dataConfigOps.planeDataConfig);
				break;

			case 4:
			{
				sweet::SphereData_Physical u(dataConfigOps.sphereDataConfig);
				sweet::SphereData_Physical v(dataConfigOps.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				dataConfigOps.ops.vrtdiv_to_uv(dataConfigOps.prog.vrt, dataConfigOps.prog.div, u, v);
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(u, dataConfigOps.planeDataConfig);
				break;
			}

			case 5:
			{
				sweet::SphereData_Physical u(dataConfigOps.prog.vrt.sphereDataConfig);
				sweet::SphereData_Physical v(dataConfigOps.prog.vrt.sphereDataConfig);

				// Don't use Robert, since we're not interested in the Robert formulation here
				dataConfigOps.ops.vrtdiv_to_uv(dataConfigOps.prog.vrt, dataConfigOps.prog.div, u, v);
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(v, dataConfigOps.planeDataConfig);
				break;
			}

			case 6:
			case 7:
			case 8:
			{
				sweet::SphereData_Spectral anal_solution_phi_pert(dataConfigOps.sphereDataConfig);
				sweet::SphereData_Spectral anal_solution_vrt(dataConfigOps.sphereDataConfig);
				sweet::SphereData_Spectral anal_solution_div(dataConfigOps.sphereDataConfig);

				sphereBenchmarks.benchmark->getInitialState(
						anal_solution_phi_pert,
						anal_solution_vrt,
						anal_solution_div
					);

				switch (id)
				{
				case 6:
					dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(dataConfigOps.prog.phi_pert - anal_solution_phi_pert, dataConfigOps.planeDataConfig);
					break;

				case 7:
					dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(dataConfigOps.prog.vrt - anal_solution_vrt, dataConfigOps.planeDataConfig);
					break;

				case 8:
					dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(dataConfigOps.prog.div - anal_solution_div, dataConfigOps.planeDataConfig);
					break;
				}
			}
		}

		double viz_min = dataConfigOps.vis_plane_data.physical_reduce_min();
		double viz_max = dataConfigOps.vis_plane_data.physical_reduce_max();

		viz_max = std::max(std::abs(viz_max), std::abs(viz_min));
		viz_min = -viz_max;

		*o_viz_min = viz_min;
		*o_viz_max = viz_max;


		*o_dataArray = &dataConfigOps.vis_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		std::ostringstream ss;
		std::string sep = ",  ";

		o_replace_commas_with_newline = true;

		if (shackParallelization->useMPI) {
			ss << "Rank=" << shackParallelization->mpiRank;
			ss << ",";
		}

		const char *fields_array[] = {
				"phi_pert",
				"vrt",
				"div",
				"h",
				"u",
				"v",
				"phi diff t0",
				"u diff t0",
				"v diff t0",
		};

		int days = ((int)(shackTimestepControl->current_simulation_time/(60.0*60.0*24.0)));
		int hours = ((int)(shackTimestepControl->current_simulation_time/(60.0*60.0)))%24;
		int minutes = ((int)(shackTimestepControl->current_simulation_time/(60.0)))%60;
		int seconds = ((int)shackTimestepControl->current_simulation_time) % 60;

		ss << "time="
				<< days
				<< "d-"
				<< (hours < 10 ? "0" : "") << hours
				<< "h:"
				<< (minutes < 10 ? "0" : "") << minutes
				<< "m:"
				<< (seconds < 10 ? "0" : "") << seconds
				<< "s"
				<< sep;


		int id = dataConfigOps.vis_data_id % max_vis_types;
		ss << "id=" << dataConfigOps.vis_data_id << sep;
		ss << "field=" << fields_array[id] << sep;
		ss << "max=" << dataConfigOps.vis_plane_data.physical_reduce_max() << sep;
		ss << "min=" << dataConfigOps.vis_plane_data.physical_reduce_min() << sep;

		ss << "time.step.nr="	<< shackTimestepControl->current_timestep_nr << sep;
		ss << "time.step.size="	<< shackTimestepControl->current_timestep_size << sep;

		if (shackPDESWESphere->compute_diagnostics)
		{
			update_diagnostics();

			ss << "diag.total_mass="	<< diagnostics.total_mass << sep;
			ss << "diag.total_energy="	<< diagnostics.total_energy << sep;
			ss << "diag.total_potential_enstrophy="	<< diagnostics.total_potential_enstrophy;
			ss << ",  ";
		}
		ss << "FIN";

		return ss.str();
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
			dataConfigOps.vis_data_id++;
			break;

		case 'V':
			dataConfigOps.vis_data_id--;
			break;

		case 'b':
			dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id + 1) % 2;
			break;

		case 'B':
			dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id - 1) % 2;
			break;
		}
	}
#endif
};





#if 0

////////////////////////////////////////////////

/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class ProgramPDESWESphereXXX
{
public:
	sweet::SphereOperators op;
	sweet::SphereOperators op_nodealiasing;

	PDESWESphere_TimeSteppers timeSteppers;

	sweet::SphereData_Spectral prog_phi_pert;
	sweet::SphereData_Spectral prog_vrt;
	sweet::SphereData_Spectral prog_div;

#if SWEET_GUI
	sweet::PlaneData_Physical vis_plane_data;
#endif

	int vis_render_type_of_primitive_id = 1;

	// was the output of the time step already done for this simulation state?
	double timestep_nr_last_output_simtime;


	void reset()
	{
		StopwatchBox::getInstance().main_setup.start();

		shackDict.reset();
		shackIOData->output_time_scale = 1.0/(60.0*60.0);

		shackIOData->output_next_sim_seconds = 0;

		if (shackTimestepControl->current_timestep_size <= 0)
			SWEETError("Only fixed time step size supported");

		// use dealiased physical space for setup
		sphereBenchmarks.setup(shackDict, op);
		sphereBenchmarks.timestepper->getInitialState(prog_phi_pert, prog_vrt, prog_div);

		/*
		 * SETUP time steppers
		 */
		timeSteppers.setup(shackDict.disc.timestepping_method,
				op, shackDict);

		std::cout << "[MULE] timestepper_string_id: " << timeSteppers.timestepper->string_id() << std::endl;

		update_diagnostics();

		diagnostics.backup_reference();

		StopwatchBox::getInstance().main_setup.stop();

		// start at one second in the past to ensure output at t=0
		timestep_nr_last_output_simtime = shackTimestepControl->current_simulation_time-1.0;

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
		if (shackIOData->output_each_sim_seconds >= 0)
			_timestepDoOutput();
	}



public:
	bool should_quit()
	{
		if (shackTimestepControl->max_timesteps_nr != -1 && shackTimestepControl->max_timesteps_nr <= shackTimestepControl->current_timestep_nr)
			return true;

		double diff = std::abs(shackTimestepControl->max_simulation_time - shackTimestepControl->current_simulation_time);

		if (	shackTimestepControl->max_simulation_time != -1 &&
				(
						shackTimestepControl->max_simulation_time <= shackTimestepControl->current_simulation_time	||
						diff/shackTimestepControl->max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
				)
			)
			return true;

		return false;
	}



	void runTimestep()
	{
#if SWEET_GUI
		if (shackDict.misc.gui_enabled && shackDict.misc.normal_mode_analysis_generation == 0)
			timestepHandleOutput();
#endif

		if (shackTimestepControl->current_simulation_time + shackTimestepControl->current_timestep_size > shackTimestepControl->max_simulation_time)
			shackTimestepControl->current_timestep_size = shackTimestepControl->max_simulation_time - shackTimestepControl->current_simulation_time;

		timeSteppers.timestepper->runTimestep(
				prog_phi_pert, prog_vrt, prog_div,
				shackTimestepControl->current_timestep_size,
				shackTimestepControl->current_simulation_time
			);



		// Apply viscosity at posteriori, for all methods explicit diffusion for non spectral schemes and implicit for spectral

		if (shackPDESWESphere->viscosity != 0 && shackDict.misc.use_nonlinear_only_visc == 0)
		{
			///prog_vrt = op.implicit_diffusion(prog_vrt, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->sphere_radius);
			///prog_div = op.implicit_diffusion(prog_div, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->sphere_radius);
			///prog_phi_pert = op.implicit_diffusion(prog_phi_pert, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->sphere_radius);
			prog_vrt = op.implicit_hyperdiffusion(prog_vrt, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->viscosity_order, shackPDESWESphere->sphere_radius);
			prog_div = op.implicit_hyperdiffusion(prog_div, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->viscosity_order, shackPDESWESphere->sphere_radius);
			prog_phi_pert = op.implicit_hyperdiffusion(prog_phi_pert, shackTimestepControl->current_timestep_size*shackPDESWESphere->viscosity, shackPDESWESphere->viscosity_order, shackPDESWESphere->sphere_radius);
		}


		// advance time step and provide information to parameters
		shackTimestepControl->current_simulation_time += shackTimestepControl->current_timestep_size;
		shackTimestepControl->current_timestep_nr++;

#if SWEET_GUI
		timestepHandleOutput();
#endif
	}
};

#endif


#if 0

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

	if (shackIOData->verbosity > 3)
		std::cout << " + setup SH sphere transformations..." << std::endl;

	sphereDataConfigInstance.setupAuto(
			shackDict.disc.space_res_physical,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackIOData->verbosity,
			shackDict.parallelization.num_threads_space
		);

	int res_physical_nodealias[2] = {
			2*shackDict.disc.space_res_spectral[0],
			shackDict.disc.space_res_spectral[1]
		};

	if (shackIOData->verbosity > 3)
		std::cout << " + setup SH sphere transformations (nodealiasing)..." << std::endl;

	sphereDataConfigInstance_nodealiasing.setupAuto(
			res_physical_nodealias,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackIOData->verbosity,
			shackDict.parallelization.num_threads_space
		);


#if SWEET_GUI
	if (shackIOData->verbosity > 3)
		std::cout << " + setup FFT plane transformations..." << std::endl;

	planeDataConfigInstance.setupAutoSpectralSpaceFromPhysical(shackDict.disc.space_res_physical, shackDict.misc.reuse_spectral_transformation_plans);
#endif

	std::ostringstream buf;
	buf << std::setprecision(14);

	if (shackIOData->verbosity > 3)
		std::cout << " + setup finished" << std::endl;

#if SWEET_MPI
	std::cout << "Hello from MPI rank: " << mpi_rank << std::endl;

	// only start simulation and time stepping for first rank
	if (mpi_rank > 0)
	{
		/*
		 * Deactivate all output for ranks larger than the current one
		 */
		shackIOData->verbosity = 0;
	#if !SWEET_XBRAID
		shackIOData->output_each_sim_seconds = -1;
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

			shackIOData->output_time_scale = 1.0/(60.0*60.0);

			//SphereOperators op(sphereDataConfig, shackPDESWESphere->plane_domain_size, shackDict.disc.space_use_spectral_basis_diffs);
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
					frac = shackTimestepControl->current_timestep_size / shackDict.parareal.coarse_timestep_size;
				else
					frac = shackTimestepControl->current_timestep_size / (shackTimestepControl->max_simulation_time / shackDict.parareal.coarse_slices );
				for (int j = 0; j < 2; j++)
					N_spectral[j] = std::max(4, int(shackDict.disc.space_res_spectral[j] * frac));

				sphereDataConfigs.push_back(new sweet::SphereData_Config);
				sphereDataConfigs.back()->setupAuto(
						N_physical,
						N_spectral,
						shackDict.misc.reuse_spectral_transformation_plans,
						shackIOData->verbosity,
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


			PDESWESphere_TimeSteppers* timeSteppersFine = new PDESWESphere_TimeSteppers;
			PDESWESphere_TimeSteppers* timeSteppersCoarse = new PDESWESphere_TimeSteppers;

			/*
			 * Allocate parareal controller and provide class
			 * which implement the parareal features
			 */
			Parareal_Controller<PDESWESphere_TimeSteppers, 3> parareal_Controller(	&shackDict,
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

			shackIOData->output_time_scale = 1.0/(60.0*60.0);

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
							shackIOData->verbosity,
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
			int nt = (int) (shackTimestepControl->max_simulation_time / shackTimestepControl->current_timestep_size);
                        if (nt * shackTimestepControl->current_timestep_size < shackTimestepControl->max_simulation_time - 1e-10)
				nt++;
			///sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackTimestepControl->max_simulation_time, nt, &shackDict, sphereDataConfig, &op);
			sweet_BraidApp app(MPI_COMM_WORLD, mpi_rank, 0., shackTimestepControl->max_simulation_time, nt, &shackDict, sphereDataConfigs, ops);


			if ( shackDict.xbraid.xbraid_run_wrapper_tests)
			{
				app.setup();

				BraidUtil braid_util;
				int test = braid_util.TestAll(&app, comm, stdout, 0., shackTimestepControl->current_timestep_size, shackTimestepControl->current_timestep_size * 2);
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
				simulationSWE->timestepHandleOutput();
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
					simulationSWE->timestepHandleOutput();

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

				if (shackIOData->verbosity > 0)
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
				simulationSWE->timestepHandleOutput();
			}

#if SWEET_MPI
			// Start counting time
			if (mpi_rank == 0)
#endif
			{
				if (shackIOData->output_file_name.size() > 0)
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
		std::cout << "[MULE] shackTimestepControl->current_timestep_nr: " << shackTimestepControl->current_timestep_nr << std::endl;
		std::cout << "[MULE] shackTimestepControl->current_timestep_size: " << shackTimestepControl->current_timestep_size << std::endl;
		std::cout << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "* Other timing information (derived)" << std::endl;
		std::cout << "***************************************************" << std::endl;
		std::cout << "[MULE] simulation_benchmark_timings.time_per_time_step (secs/ts): " << StopwatchBox::getInstance().main_timestepping()/(double)shackTimestepControl->current_timestep_nr << std::endl;
	}

#if SWEET_MPI
	MPI_Finalize();
#endif

	return 0;
}
#endif

#endif
