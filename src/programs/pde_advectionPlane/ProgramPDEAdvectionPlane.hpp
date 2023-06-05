/*
 * 		Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_ADVECTIONPLANE_PROGRAMPDEADVECTIONPLANE_HPP_
#define SRC_PROGRAMS_PDE_ADVECTIONPLANE_PROGRAMPDEADVECTIONPLANE_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include "PDEAdvectionPlaneTimeSteppers.hpp"
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/plane/Plane.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>

// Benchmarks
#include "PDEAdvectionPlaneBenchmarksCombined.hpp"

// Time steppers

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
#endif


class ProgramPDEAdvectionPlane
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

		sweet::PlaneData_Config planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral prog_h;
		sweet::PlaneData_Spectral prog_h_t0;	// at t0
		sweet::PlaneData_Spectral prog_u;
		sweet::PlaneData_Spectral prog_v;


		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			planeDataConfig.setupAuto(*i_shackPlaneDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeDataConfig);

			ops.setup(planeDataConfig, *i_shackPlaneDataOps);
			ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(ops);

			prog_h.setup(planeDataConfig);
			prog_h_t0.setup(planeDataConfig);
			prog_u.setup(planeDataConfig);
			prog_v.setup(planeDataConfig);

			return true;
		}

		void clear()
		{
			prog_h.clear();
			prog_h_t0.clear();
			prog_u.clear();
			prog_v.clear();

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// time integrators
	PDEAdvectionPlaneTimeSteppers timeSteppers;

	// Handler to all benchmarks
	PDEAdvectionPlaneBenchmarksCombined planeBenchmarksCombined;

	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;
	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDEAdvectionPlaneTimeDiscretization *shackTimeDisc;


#if SWEET_GUI
	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical vis_plane_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 0;

	// Which primitive to use for rendering
	int vis_data_id = 0;

#endif

public:
	ProgramPDEAdvectionPlane(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackPlaneDataOps(nullptr),
		shackIOData(nullptr),
		shackTimestepControl(nullptr),
		shackTimeDisc(nullptr)
	{
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN(shackProgArgDict);
	}



	bool setup_1_shackRegistration()
	{
		/*
		 * SHACK: Register classes which we require
		 */
		shackPlaneDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		shackTimestepControl = shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		shackTimeDisc = shackProgArgDict.getAutoRegistration<ShackPDEAdvectionPlaneTimeDiscretization>();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register other things before parsing program arguments
		 */
		planeBenchmarksCombined.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);

		timeSteppers.shackRegistration(shackProgArgDict);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

		shackProgArgDict.processHelpArguments();
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_1_shackRegistration()
	{
		shackPlaneDataOps = nullptr;
		shackTimestepControl = nullptr;
		shackIOData = nullptr;
		shackTimeDisc = nullptr;

		planeBenchmarksCombined.clear();
		timeSteppers.clear();
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

	bool setup_3_dataOpsEtc()
	{
		/*
		 * Setup Plane Data Config & Operators
		 */
		dataConfigOps.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(dataConfigOps);

		/*
		 * After we setup the plane, we can setup the time steppers and their buffers
		 */
		timeSteppers.setup(shackProgArgDict, dataConfigOps.ops);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(timeSteppers);

#if SWEET_GUI
		vis_plane_data.setup(dataConfigOps.planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackProgArgDict.printShackData();

		planeBenchmarksCombined.setupInitialConditions(
				dataConfigOps.prog_h,
				dataConfigOps.prog_u,
				dataConfigOps.prog_v,
				dataConfigOps.ops,
				shackProgArgDict
			);
		ERROR_CHECK_WITH_FORWARD_AND_COND_RETURN_BOOLEAN(planeBenchmarksCombined);

		dataConfigOps.prog_h_t0 = dataConfigOps.prog_h;

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
		vis_plane_data.clear();
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

	void printSimulationErrors()
	{
		std::cout << "Error compared to initial condition" << std::endl;
		std::cout << "Lmax error: " << (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toPhys().physical_reduce_max_abs() << std::endl;
		std::cout << "RMS error: " << (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toPhys().physical_reduce_rms() << std::endl;
	}

	double getErrorLMaxOnH()
	{
		return (dataConfigOps.prog_h_t0-dataConfigOps.prog_h).toPhys().physical_reduce_max_abs();
	}

	virtual ~ProgramPDEAdvectionPlane()
	{
		clear();
	}

	bool runTimestep()
	{
		shackTimestepControl->timestepHelperStart();

		timeSteppers.master->runTimestep(
				dataConfigOps.prog_h, dataConfigOps.prog_u, dataConfigOps.prog_v,
				shackTimestepControl->current_timestepSize,
				shackTimestepControl->current_simulation_time
			);

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
	void vis_post_frame_processing(
		int i_num_iterations
	)
	{
		if (shackTimestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				runTimestep();
	}

	void vis_getDataArray(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *vis_reset
	)
	{
		*o_render_primitive_id = vis_render_type_of_primitive_id;

		int id = vis_data_id % 3;
		switch (id)
		{
		case 0:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(dataConfigOps.prog_h, vis_plane_data);
			break;

		case 1:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(dataConfigOps.prog_u, vis_plane_data);
			break;

		case 2:
			sweet::Convert_PlaneDataSpectral_To_PlaneDataPhysical::convert(dataConfigOps.prog_v, vis_plane_data);
			break;
		}

		*o_dataArray = &vis_plane_data;
		*o_aspect_ratio = 1;
	}


	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		const char* description = "";
		int id = vis_data_id % 3;

		switch (id)
		{
		default:
		case 0:
			description = "H";
			break;

		case 1:
			description = "u";
			break;

		case 2:
			description = "v";
			break;
		}

		static char title_string[2048];

		sprintf(title_string,
#if SWEET_MPI
				"Rank %i - "
#endif
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				-1,	// TODO: mpi_rank,
#endif
				shackTimestepControl->current_simulation_time,
				shackTimestepControl->current_simulation_time/(60.0*60.0*24.0),
				shackTimestepControl->current_timestep_nr,
				shackTimestepControl->current_timestepSize,
				description,
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

		case 'b':
			vis_render_type_of_primitive_id = (vis_render_type_of_primitive_id + 1) % 2;
			break;
		}
	}
#endif
};




#endif
