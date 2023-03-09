/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 * 
 * MULE_SCONS_OPTIONS: --plane-spectral-space=enable
 */

#include <ostream>
#include <sstream>

// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space required"
#endif

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Parse program arguments
#include <sweet/core/ProgramArguments.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/plane/Plane.hpp>

// Our shack directory to store different objects and get them back later on
#include <sweet/core/shacks/ShackDictionary.hpp>

#include <sweet/core/Stopwatch.hpp>

#include <sweet/gui/VisSweet.hpp>



class ProgramPlaneSpectralVisualization
#if SWEET_GUI
		:	public SimulationGUICallbacks
#endif
{
public:
	sweet::ErrorBase error;
	sweet::ProgramArguments programArguments;

	class Data
	{
	public:
		sweet::ErrorBase error;

		sweet::PlaneData_Config planeDataConfig;
		sweet::PlaneOperators ops;

		sweet::PlaneData_Spectral tmp;
		sweet::PlaneData_Physical tmp_phys;

		bool setup(sweet::ShackPlaneDataOps *i_shackPlaneDataOps)
		{
			/*
			 * Setup Plane Data Config & Operators
			 */
			if (!planeDataConfig.setupAuto(*i_shackPlaneDataOps))
				return error.forwardWithPositiveReturn(planeDataConfig.error);

			if (!ops.setup(
					planeDataConfig,
					*i_shackPlaneDataOps
				))
				return error.forwardWithPositiveReturn(ops.error);

			tmp.setup(planeDataConfig);
			tmp_phys.setup(planeDataConfig);

			return true;
		}

		void clear()
		{
			tmp.clear();
			tmp_phys.clear();

			ops.clear();
			planeDataConfig.clear();
		}
	};

	// Data
	Data data;
	
	sweet::ShackDictionary shackDict;
	sweet::ShackPlaneDataOps *shackPlaneDataOps;


	// Data to visualize is stored to this variable
	sweet::PlaneData_Physical vis_plane_data;

	// Which primitive to use for rendering
	int vis_dataId = 0;

	std::string vis_description;


	Stopwatch stopwatch;

	/*
	 * Setup main
	 */
	int prog_argc;
	char *const * const prog_argv;

public:
	ProgramPlaneSpectralVisualization(
			int i_argc,
			char *const * const i_argv
	)	:
		prog_argc(i_argc),
		prog_argv(i_argv)
	{
	}


	/*
	 * Clear all data
	 */
	void clear()
	{
		programArguments.clear();
		shackDict.clear();

		data.clear();

#if SWEET_GUI
		vis_plane_data.clear();
#endif
	}

	/*
	 * Chekc if help should be printed a do so
	 */
	void checkAndPrintHelp()
	{
		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
		}
	}


	bool setup()
	{
		/*
		 * Parse program arguments
		 */
		programArguments.setup(prog_argc, prog_argv);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(programArguments);

		/*
		 * SHACK: Register classes which we require
		 */
		shackPlaneDataOps = shackDict.getAutoRegistration<sweet::ShackPlaneDataOps>();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackDict);

		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			shackDict.printProgramArguments();
			return false;
		}

		/*
		 * SHACK: Process arguments
		 */
		shackDict.processProgramArguments(programArguments);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(shackDict);

		/*
		 * Setup Plane Data Config & Operators
		 */
		data.setup(shackPlaneDataOps);
		ERROR_CHECK_WITH_RETURN_BOOLEAN(data);


#if SWEET_GUI
		vis_plane_data.setup(data.planeDataConfig);
#endif

		std::cout << "Printing shack information:" << std::endl;
		shackDict.printShackData();

		/*
		 * Finish registration & getting class interfaces so that nobody can do some
		 * strange things with this anymore
		 */
		shackDict.closeRegistration();
		shackDict.closeGet();

		/*
		 * Now we should check that all program arguments have really been parsed
		 */
		programArguments.checkAllArgumentsProcessed();
		ERROR_CHECK_WITH_RETURN_BOOLEAN(programArguments);

		return true;
	}

	bool reset()
	{
		clear();
		setup();

		return !error.exists();
	}


	bool runTimestep()
	{
		return true;
	}


public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
	}



public:
	bool should_quit()
	{
		return false;
	}


	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
		int i_num_iterations
	)
	{
	}



	void vis_get_vis_data_array(
			const sweet::PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_vis_min,
			double *o_vis_max,
			bool *vis_reset
	)
	{
		data.tmp.spectral_set_zero();
		int spec_array[][4] =
		{
				// j, i, re, im
				{0, 0, 1, 0},
				{0, 0, 0, 1},
				{0, 0, 1, 1},

				{1, 0, 1, 0},
				{1, 0, 0, 1},
				{1, 0, 1, 1},

				{2, 0, 1, 0},
				{2, 0, 0, 1},
				{2, 0, 1, 1},

				{0, 1, 1, 0},
				{0, 1, 0, 1},
				{0, 1, 1, 1},

				{1, 1, 1, 0},
				{1, 1, 0, 1},
				{1, 1, 1, 1},

				{(int)shackPlaneDataOps->space_res_physical[1]-1, 0, 1, 0},
				{(int)shackPlaneDataOps->space_res_physical[1]-1, 0, 0, 1},
				{(int)shackPlaneDataOps->space_res_physical[1]-1, 0, 1, 1},

				{(int)shackPlaneDataOps->space_res_physical[1]-2, 0, 1, 0},
				{(int)shackPlaneDataOps->space_res_physical[1]-2, 0, 0, 1},
				{(int)shackPlaneDataOps->space_res_physical[1]-2, 0, 1, 1},

				{0, (int)shackPlaneDataOps->space_res_physical[0]/2, 1, 0},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2, 0, 1},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2, 1, 1},

				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-1, 1, 0},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-1, 0, 1},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-1, 1, 1},

				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-2, 1, 0},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-2, 0, 1},
				{0, (int)shackPlaneDataOps->space_res_physical[0]/2-2, 1, 1},
		};

		int id = vis_dataId % (sizeof(spec_array)/sizeof(spec_array[0]));

		std::ostringstream ss;
		ss << "spec_coord (j, i) = (" << spec_array[id][0] << ", " << spec_array[id][1] << ")";
		ss << ", ";
		ss << "value = (" << spec_array[id][2] << ", " << spec_array[id][3] << "i)";
		vis_description = ss.str();

		data.tmp.spectral_set(
				spec_array[id][0],
				spec_array[id][1],
				{
						(double)spec_array[id][2],
						(double)spec_array[id][3]
				}
			);

		data.tmp_phys = data.tmp.toPhys();

		*o_dataArray = &data.tmp_phys;

		*o_aspect_ratio = shackPlaneDataOps->plane_domain_size[1] / shackPlaneDataOps->plane_domain_size[0];
	}


	/**
	 * return status string for window title
	 */
	const std::string vis_getStatusString(bool &o_replace_commas_with_newline)
	{
		o_replace_commas_with_newline = false;
		return vis_description;
	}



	void vis_pause()
	{
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			vis_dataId++;
			break;

		case 'V':
			vis_dataId--;
			break;
		}
	}


};


int main(int i_argc, char *i_argv[])
{
	ProgramPlaneSpectralVisualization simulation(i_argc, i_argv);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	simulation.setup();
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);


	VisSweet visSweet(simulation);
	ERROR_CHECK_WITH_PRINT_AND_RETURN_EXIT(simulation);

	std::cout << "FIN" << std::endl;
	return 0;
}
