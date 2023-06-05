/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTIONSPHERE_PROGRAMPDEADVECTIONSPHERE_HPP_
#define SRC_PROGRAMS_PDE_ADVECTIONSPHERE_PROGRAMPDEADVECTIONSPHERE_HPP_


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
#include "benchmarks/ShackPDEAdvectionSphereBenchmarks.hpp"
#include "ShackPDEAdvectionSphere.hpp"

// Benchmarks
#include "PDEAdvectionSphereBenchmarksCombined.hpp"

// Time steppers
#include "PDEAdvectionSphereTimeSteppers.hpp"

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
	#include <sweet/core/plane/Plane.hpp>
	#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#endif

#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/sphere/SphereData_DebugContainer.hpp>



class ProgramPDEAdvectionSphere
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

		std::vector<sweet::SphereData_Spectral> prog_vec;
		sweet::SphereData_Physical vel_u;
		sweet::SphereData_Physical vel_v;

		std::vector<sweet::SphereData_Spectral> prog_vec_t0;

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
				int i_num_vec_elements
		)
		{
			/*
			 * Setup Sphere Data Config & Operators
			 */
			sphereDataConfig.setupAuto(i_shackSphereDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(sphereDataConfig);

			ops.setup(&sphereDataConfig, i_shackSphereDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			vel_u.setup(sphereDataConfig);
			vel_v.setup(sphereDataConfig);

			prog_vec.resize(i_num_vec_elements);
			prog_vec_t0.resize(i_num_vec_elements);
			for (int i = 0; i < i_num_vec_elements; i++)
			{
				prog_vec[i].setup(sphereDataConfig);
				prog_vec_t0[i].setup(sphereDataConfig);
			}

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
			for (std::size_t i = 0; i < prog_vec.size(); i++)
			{
				prog_vec[i].clear();
				prog_vec_t0[i].clear();
			}

			prog_vec.clear();
			prog_vec_t0.clear();

			vel_u.clear();
			vel_v.clear();

			ops.clear();
			sphereDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// time integrators
	PDEAdvectionSphereTimeSteppers timeSteppers;

	// Handler to all benchmarks
	PDEAdvectionSphereBenchmarksCombined benchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDEAdvectionSphereTimeDiscretization *shackTimeDisc;
	ShackPDEAdvectionSphereBenchmarks *shackBenchmarks;
	ShackPDEAdvectionSphere *shackPDEAdvectionSphere;



public:
	ProgramPDEAdvectionSphere(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr),
		shackBenchmarks(nullptr),
		shackPDEAdvectionSphere(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
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
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphereTimeDiscretization>();
		shackBenchmarks = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphereBenchmarks>();
		shackPDEAdvectionSphere = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionSphere>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */

		/*
		 * Setup benchmarks
		 */
		benchmarksCombined.setup_1_registerAllBenchmark();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);

		benchmarksCombined.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Setup time steppers
		 */
		timeSteppers.setup_1_registerAllTimesteppers();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		timeSteppers.setup_2_shackRegistration(&shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		/*
		 * Process HELP arguments
		 */
		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

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

	bool setup_3_dataOpsEtc()
	{
		/*
		 * BENCHMARK: Detect particular benchmark to use
		 */
		benchmarksCombined.setup_3_benchmarkDetection();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);

		/*
		 * Setup benchmark itself
		 */
		benchmarksCombined.setup_4_benchmarkSetup_1_withoutOps();

		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphereDataOps, benchmarksCombined.benchmark->getNumPrognosticFields());
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

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
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);


		/*
		 * Load initial state of benchmark
		 */
		benchmarksCombined.benchmark->getInitialState(
				dataConfigOps.prog_vec,
				dataConfigOps.vel_u,
				dataConfigOps.vel_v
			);

		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(benchmarksCombined);


		/*
		 * Backup data at t=0
		 */
		dataConfigOps.prog_vec_t0 = dataConfigOps.prog_vec;

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

		if (!setup_3_dataOpsEtc())
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

	double getErrorLMaxOnH()
	{
		return (dataConfigOps.prog_vec_t0[0]-dataConfigOps.prog_vec[0]).toPhys().physical_reduce_max_abs();
	}

	double getErrorRMSOnH()
	{
		return (dataConfigOps.prog_vec_t0[0]-dataConfigOps.prog_vec[0]).toPhys().physical_reduce_rms();
	}

	void printSimulationErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << getErrorLMaxOnH() << std::endl;
		std::cout << "RMS error: " << getErrorRMSOnH() << std::endl;
	}

	virtual ~ProgramPDEAdvectionSphere()
	{
		clear();
	}


	bool runTimestep()
	{
		sweet::SphereData_DebugContainer::clear();

		shackTimestepControl->timestepHelperStart();


		timeSteppers.timestepper->runTimestep(
				dataConfigOps.prog_vec, dataConfigOps.vel_u, dataConfigOps.vel_v,
				shackTimestepControl->current_timestepSize,
				shackTimestepControl->current_simulation_time
			);

		if (shackBenchmarks->getVelocities)
		{
			/*
			 * Update velocities just for sake of the correction visualization
			 */
			shackBenchmarks->getVelocities(
					dataConfigOps.vel_u,
					dataConfigOps.vel_v,
					shackTimestepControl->current_simulation_time + shackTimestepControl->current_timestepSize,
					shackBenchmarks->getVelocitiesUserData
				);
		}

		shackTimestepControl->timestepHelperEnd();

		if (shackIOData->verbosity > 10)
		{
			std::cout << "ts_nr=" << shackTimestepControl->current_timestep_nr << ", t=" << shackTimestepControl->current_simulation_time*shackIOData->output_time_scale_inv << std::endl;
			std::cout << "error:" << getErrorLMaxOnH() << std::endl;
		}

		return true;
	}



	bool should_quit()
	{
		return shackTimestepControl->isFinalTimestepReached();
	}



#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				runTimestep();
	}


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
		*o_vis_render_type_of_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
		*o_bogus_data = &dataConfigOps.sphereDataConfig;

		if (dataConfigOps.vis_data_id < 0)
		{
			int id = -dataConfigOps.vis_data_id-1;

#if 0
			if (id <  (int)SphereData_DebugContainer().size())
			{
				sweet::SphereData_DebugContainer::DataContainer &d = SphereData_DebugContainer().container_data()[id];
				if (d.is_spectral)
					vis_plane_data = Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(d.data_spectral, sphereDataConfig);
				else
					vis_plane_data = Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(d.data_physical, sphereDataConfig);

				*o_dataArray = &vis_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
#else
			if (id <  (int)dataConfigOps.prog_vec.size())
			{
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(
							dataConfigOps.prog_vec[id] - dataConfigOps.prog_vec_t0[id],
							dataConfigOps.planeDataConfig
						);

				*o_dataArray = &dataConfigOps.vis_plane_data;
				*o_aspect_ratio = 0.5;
				return;
			}
#endif
		}

		std::size_t id = dataConfigOps.vis_data_id;

		if (id >= 0 && id < dataConfigOps.prog_vec.size())
		{
			dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(dataConfigOps.prog_vec[id], dataConfigOps.planeDataConfig);
		}
		else if (id >= dataConfigOps.prog_vec.size() && id < dataConfigOps.prog_vec.size() + 2)
		{
			switch (id - dataConfigOps.prog_vec.size())
			{
			case 0:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(dataConfigOps.vel_u, dataConfigOps.planeDataConfig);
				break;

			case 1:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(dataConfigOps.vel_v, dataConfigOps.planeDataConfig);
				break;
			}
		}
		else if (
			dataConfigOps.prog_vec.size() == 2		&&
			id >= dataConfigOps.prog_vec.size() + 2	&&
			id < dataConfigOps.prog_vec.size() + 4
		)
		{
			sweet::SphereData_Physical u, v;
			dataConfigOps.ops.vrtdiv_to_uv(dataConfigOps.prog_vec[0], dataConfigOps.prog_vec[1], u, v);

			switch (id - dataConfigOps.prog_vec.size() - 2)
			{
			case 0:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(u, dataConfigOps.planeDataConfig);
				break;

			case 1:
				dataConfigOps.vis_plane_data = sweet::Convert_SphereDataPhysical_To_PlaneDataPhysical::physical_convert(v, dataConfigOps.planeDataConfig);
				break;
			}
		}
		else
		{
			SWEETDebugAssert(dataConfigOps.vis_plane_data.physical_space_data != nullptr);
			SWEETDebugAssert(dataConfigOps.vis_plane_data.planeDataConfig != nullptr);
			dataConfigOps.vis_plane_data.physical_set_zero();
		}

		*o_dataArray = &dataConfigOps.vis_plane_data;
		*o_aspect_ratio = 0.5;
	}



	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		std::string description = "";

		bool found = false;
		if (dataConfigOps.vis_data_id < 0)
		{
			int id = -dataConfigOps.vis_data_id-1;

#if 0
			if (id < (int)SphereData_DebugContainer().size())
			{
				description = std::string("DEBUG_")+SphereData_DebugContainer().container_data()[id].description;
				found = true;
			}
#else
			if (id <  (int)dataConfigOps.prog_vec.size())
			{
				std::ostringstream msg;
				msg << "DIFF prog. field " << id;
				description = msg.str();
			}
#endif
		}

		if (!found)
		{
			std::size_t id = dataConfigOps.vis_data_id;

			if (id >= 0 && id < dataConfigOps.prog_vec.size())
			{
				std::ostringstream msg;
				msg << "Prog. field " << id;
				description = msg.str();
			}
			else if (id >= dataConfigOps.prog_vec.size() && id < dataConfigOps.prog_vec.size() + 2)
			{
				switch (id - dataConfigOps.prog_vec.size())
				{
				case 0:
					description = "u velocity";
					break;

				case 1:
					description = "v velocity";
					break;
				}
			}
			else if (
					dataConfigOps.prog_vec.size() == 2		&&
					id >= dataConfigOps.prog_vec.size() + 2	&&
					id < dataConfigOps.prog_vec.size() + 4
			)
			{
				switch (id - dataConfigOps.prog_vec.size() - 2)
				{
				case 0:
					description = "prognostic field: u velocity";
					break;

				case 1:
					description = "prognostic field: v velocity";
					break;
				}
			}
			else
			{
				description = "field doesn't exist";
			}
		}

		static char title_string[2048];

		//sprintf(title_string, "Time (days): %f (%.2f d), Timestep: %i, timestep size: %.14e, Vis: %s, Mass: %.14e, Energy: %.14e, Potential Entrophy: %.14e",
		sprintf(title_string,
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e "
				","
				"Colorscale: lowest [Blue... green ... red] highest",
				shackTimestepControl->current_simulation_time,
				shackTimestepControl->current_simulation_time/(60.0*60.0*24.0),
				shackTimestepControl->current_timestep_nr,
				shackTimestepControl->current_timestepSize,
				description.c_str(),
				dataConfigOps.vis_plane_data.physical_reduce_max(),
				dataConfigOps.vis_plane_data.physical_reduce_min()
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
			dataConfigOps.vis_data_id++;
			break;

		case 'V':
			dataConfigOps.vis_data_id--;
			break;

		case 'b':
		case 'B':
			dataConfigOps.vis_render_type_of_primitive_id = (dataConfigOps.vis_render_type_of_primitive_id + 1) % 2;
			break;
		}
	}
#endif
};



#endif /* SRC_PROGRAMS_PDE_ADVECTIONSPHERE_PROGRAMPDEADVECTIONSPHERE_HPP_ */
