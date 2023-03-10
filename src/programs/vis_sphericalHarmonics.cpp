/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 * MULE_SCONS_OPTIONS: --gui=enable
 */


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Include everything we need for simulations on the sphere and plane
#include <sweet/core/sphere/Sphere.hpp>
#include <sweet/core/plane/Plane.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>

#if SWEET_GUI
	#include <sweet/gui/VisSweet.hpp>
	#include <sweet/core/plane/Plane.hpp>
	#include <sweet/core/shacksShared/ShackPlaneDataOps.hpp>
#endif

#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include <sweet/core/sphere/SphereData_DebugContainer.hpp>

#include "pde_sweSphere/PDESWESphereBenchmarksCombined.hpp"

#include <sweet/core/SWEETError.hpp>
#include <sweet/core/sphere/SphereData_DebugContainer.hpp>


/*
 * This allows running REXI including Coriolis-related terms but just by setting f to 0
 */


class ProgramVisSphericalHarmonics
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

		sweet::SphereData_Spectral sphereDataForVisualization;


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

			sphereDataForVisualization.setup(sphereDataConfig);

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
			sphereDataForVisualization.clear();

			ops.clear();
			sphereDataConfig.clear();
		}
	};

	// Simulation data
	DataConfigOps dataConfigOps;

	// Diagnostics measures
	int last_timestep_nr_update_diagnostics = -1;


#if SWEET_GUI && 0
	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical vis_plane_data;

	// Which primitive to use for rendering
	int vis_render_type_of_primitive_id = 1;

	// Which primitive to use for rendering
	int vis_data_id = 0;
#endif


	// was the output of the time step already done for this simulation state?
	double timestep_last_output_simtime;

	int mode_m = 0;
	int mode_n = 0;

	bool vis_reset = false;


	/*
	 * Shack directory and shacks to work with
	 */
	sweet::ShackProgArgDictionary shackProgArgDict;
	sweet::ShackSphereDataOps *shackSphereDataOps;
	sweet::ShackIOData *shackIOData;



public:
	ProgramVisSphericalHarmonics(
			int i_argc,
			char *const * const i_argv
	)	:
		shackProgArgDict(i_argc, i_argv),
		shackSphereDataOps(nullptr),
		shackIOData(nullptr)
	{
	}

	~ProgramVisSphericalHarmonics()
	{
		clear();
	}


	bool setup_1_shackRegistration()
	{
		/*
		 * Setup argument parsing
		 */
		shackProgArgDict.setup();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		/*
		 * SHACK: Register classes which we require
		 */
		shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);
		shackIOData = shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

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
		shackIOData = nullptr;

		shackProgArgDict.clear();
	}

	bool setup_2_processArguments()
	{
		/*
		 * SHACK: Process arguments
		 */
		shackProgArgDict.processProgramArguments();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackProgArgDict);

		return true;
	}

	void clear_2_processArguments()
	{
		shackProgArgDict.clear();
	}

	bool setup_3_data()
	{
		/*
		 * Setup the data fields
		 */
		dataConfigOps.setup(shackSphereDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(dataConfigOps);

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


	void setup_mode()
	{
		if (mode_n < 0)
			mode_n = 0;

		if (mode_m < 0)
			mode_m = 0;

		if (mode_n > dataConfigOps.sphereDataConfig.spectral_modes_n_max)
			mode_n = dataConfigOps.sphereDataConfig.spectral_modes_n_max;

		if (mode_m > dataConfigOps.sphereDataConfig.spectral_modes_m_max)
			mode_m = dataConfigOps.sphereDataConfig.spectral_modes_m_max;

		if (mode_m > mode_n)
			mode_m = mode_n;

		dataConfigOps.sphereDataForVisualization.spectral_set_zero();
		std::complex<double> val = 1;
		dataConfigOps.sphereDataForVisualization.spectral_set(mode_n, mode_m, val);
	}

public:
	bool should_quit()
	{
		return false;
	}


	bool runTimestep()
	{
		return true;
	}



#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
		runTimestep();
	}


	int max_vis_types = 9;


	void vis_getDataArray(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *o_vis_reset
	)
	{
		*o_vis_reset = vis_reset;
		vis_reset = false;

		// request rendering of sphere or plane
		*o_render_primitive_id = dataConfigOps.vis_render_type_of_primitive_id;
		*o_bogus_data = &dataConfigOps.sphereDataConfig;


		//int id = dataConfigOps.vis_data_id % max_vis_types;

		dataConfigOps.vis_plane_data = sweet::Convert_SphereDataSpectral_To_PlaneDataPhysical::physical_convert(dataConfigOps.sphereDataForVisualization, dataConfigOps.planeDataConfig);

		double vis_min = dataConfigOps.vis_plane_data.physical_reduce_min();
		double vis_max = dataConfigOps.vis_plane_data.physical_reduce_max();

		vis_max = std::max(std::abs(vis_max), std::abs(vis_min));
		vis_min = -vis_max;

		*o_vis_min = vis_min;
		*o_vis_max = vis_max;


		*o_dataArray = &dataConfigOps.vis_plane_data;
		*o_aspect_ratio = 0.5;
	}



	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		o_replace_commas_with_newline = true;
		std::ostringstream description;

		description << "Visualizing modes (N=" << mode_n << " | M=" << mode_m << ")";
		description << ",";
		description << "  Min: " << dataConfigOps.vis_plane_data.physical_reduce_max();
		description << ",";
		description << "  Max: " << dataConfigOps.vis_plane_data.physical_reduce_max();

		return description.str();
	}



	void vis_pause()
	{
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

		case 'm':
			mode_m++;
			setup_mode();
			break;

		case 'M':
			mode_m--;
			setup_mode();
			break;

		case 'n':
			mode_n++;
			setup_mode();
			break;

		case 'N':
			mode_n--;
			setup_mode();
			break;

		case 't':
			shackSphereDataOps->space_res_spectral[0] *= 2;
			shackSphereDataOps->space_res_spectral[1] *= 2;
			reset();
			break;

		case 'T':
			shackSphereDataOps->space_res_spectral[0] /= 2;
			shackSphereDataOps->space_res_spectral[1] /= 2;
			reset();
			break;
		}
	}
#endif
};



int main(int i_argc, char *i_argv[])
{
	ProgramVisSphericalHarmonics visSphericalHarmonics(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(visSphericalHarmonics);

	visSphericalHarmonics.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(visSphericalHarmonics);

#if SWEET_GUI
	if (visSphericalHarmonics.shackIOData->gui_enabled)
	{
		VisSweet visSweet(visSphericalHarmonics);
	}
	else
#endif
	{
		while (!visSphericalHarmonics.should_quit())
			visSphericalHarmonics.runTimestep();
	}

	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(visSphericalHarmonics);


	std::cout << "FIN" << std::endl;
	return 0;
}
