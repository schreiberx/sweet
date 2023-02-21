/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pdeAdvectionPlane/
 */

#include <sweet/defaultPrecompilerValues.hpp>
#include <sweet/ErrorBase.hpp>
#include <sweet/ProgramArguments.hpp>
#include <sweet/plane/Plane.hpp>

#include <sweet/shacks/ShackDictionary.hpp>
#include <sweet/shacksShared/ShackIOData.hpp>
#include <sweet/shacksShared/ShackDiagnostics.hpp>
#include <sweet/shacksShared/ShackDiscretization.hpp>
#include <sweet/shacksShared/ShackMisc.hpp>

#include "pdeAdvectionPlane/Adv_Plane_TimeSteppers.hpp"
#include "swe_plane_benchmarks/SWEPlaneBenchmarksCombined.hpp"



#if SWEET_GUI
	#include "sweet/VisSweet.hpp"
#endif

#pragma GCC warning "DEPRECATED! Try to avoid this!"
#include <sweet/SimulationVariables.hpp>
SimulationVariables _simVars;


#if 0
class PDEAdvectionPlaneParameters	:
		public sweet::ClassDictionaryInterface
{
public:
	sweet::ErrorBase error;

	/**
	 * Velocity and additional parameter for advection test cases
	 */
	double advection_velocity[2] = {0, 0};


	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--advection-x", advection_velocity[0]);
		i_pa.getArgumentValueByKey("--advection-y", advection_velocity[1]);

		return error.forwardFrom(i_pa.error);
	}

	void printProgramArguments(const std::string& i_prefix = "")
	{

	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "PDEAdvection:" << std::endl;
		std::cout << " + advection_velocity: " << advection_velocity[0] << ", " << advection_velocity[1] << std::endl;
	}

};
#endif

class SimulationInstance
{
public:
	sweet::ErrorBase error;

	sweet::ProgramArguments programArguments;

	PlaneDataConfig *planeDataConfig;

	PlaneOperators *ops;

	class SimPlaneData
	{
	public:
		PlaneData_Spectral prog_h;
		PlaneData_Spectral prog_h_t0;	// at t0
		PlaneData_Spectral prog_u;
		PlaneData_Spectral prog_v;

		SimPlaneData()
		{
		}

		void setup(PlaneDataConfig *i_planeDataConfig)
		{
			prog_h.setup(i_planeDataConfig);
			prog_h_t0.setup(i_planeDataConfig);
			prog_u.setup(i_planeDataConfig);
			prog_v.setup(i_planeDataConfig);
		}
	};

	SimPlaneData *simPlaneData;

	Adv_Plane_TimeSteppers timeSteppers;


#if SWEET_GUI
	PlaneData_Physical viz_plane_data;

	int render_primitive_id = 0;

	bool gui_active = false;
#endif

	SWEPlaneBenchmarksCombined planeBenchmarkCombined;

	double max_error_h0 = -1;


	SimulationCoefficients *sim;
	Diagnostics *diagnostics;
	Discretization *disc;
	Misc *misc;
	TimestepControl *timestepControl;
	ShackIOData *ioData;

	sweet::ClassInstanceDictionary classDict;


public:
	SimulationInstance()	:
		planeDataConfig(nullptr),
		ops(nullptr),
		sim(nullptr),
		diagnostics(nullptr),
		disc(nullptr),
		misc(nullptr),
		timestepControl(nullptr)
#if SWEET_GUI
		,viz_plane_data(planeDataConfig)
#endif
	{
	}

#if 1

	/*
	 * Setup main
	 */
	int argc;
	const char *const * argv;
	bool setup(int i_argc, char * i_argv[])
	{
		argc = i_argc;
		argv = i_argv;

		bool retval = reset();

		return retval;
	}

	bool reset()
	{
		clear();

		/*
		 * Parse program arguments
		 */
		if (!programArguments.setup(argc, argv))
		{
			error.forwardFrom(programArguments.error);
			return false;
		}

#if 1
		/*
		 * Register classes
		 */
		classDict.reset();

//		classDict.registerClassInstance<PDEAdvectionPlaneParameters>();
		classDict.registerClassInstance<SimulationCoefficients>();
		classDict.registerClassInstance<Diagnostics>();
		classDict.registerClassInstance<TimestepControl>();
		classDict.registerClassInstance<Discretization>();
		classDict.registerClassInstance<Misc>();

		if (classDict.error.exists())
		{
			error.forwardFrom(classDict.error);
			return false;
		}

		classDict.registrationOfClassInstancesFinished();


//		pdeAdvectionPlaneParameters = classDict.getClassInstance<PDEAdvectionPlaneParameters>();
		sim = classDict.getClassInstance<SimulationCoefficients>();
		diagnostics = classDict.getClassInstance<Diagnostics>();
		timestepControl = classDict.getClassInstance<TimestepControl>();
		disc = classDict.getClassInstance<Discretization>();
		misc = classDict.getClassInstance<Misc>();
		ioData = classDict.getClassInstance<ShackIOData>();

		if (classDict.error.exists())
		{
			error.forwardFrom(classDict.error);
			return false;
		}

		/*
		 * First, check for --help or -h
		 */
		if (programArguments.argumentWithKeyExists("-h") || programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			classDict.printProgramArguments();
			return false;
		}

		if (!classDict.processProgramArguments(programArguments))
		{
			error.forwardFrom(programArguments.error);
			return false;
		}

		if (!programArguments.checkAllArgumentsProcessed())
		{
			error.forwardFrom(programArguments.error);
			return false;
		}

		if (timestepControl->current_timestep_size <= 0)
		{
			error.set("Timestep size not set");
			return false;
		}

		std::cout << "Printing class content:" << std::endl;
		classDict.printClass();

		_simVars.reset();
#endif

#if 0
		planeBenchmarkCombined.setupInitialConditions(
				simPlaneData->prog_h,
				simPlaneData->prog_u,
				simPlaneData->prog_v,
				_simVars,
				*ops
			);

		simPlaneData->prog_h_t0 = simPlaneData->prog_h;

		timeSteppers.setup(
				pdePlaneDiscretization->timestepping_method,
				*ops,
				_simVars
			);

		_simVars.outputConfig();
#endif


		planeDataConfig = new PlaneDataConfig;
		planeDataConfig->setupAuto(
				disc->space_res_physical,
				disc->space_res_spectral,
				misc->reuse_spectral_transformation_plans
			);

		ops = new PlaneOperators;
		ops->setup(
				planeDataConfig,
				sim->plane_domain_size
			);

		simPlaneData = new SimPlaneData;
		simPlaneData->setup(planeDataConfig);

		return true;
	}


	void clear()
	{
		delete simPlaneData;
		delete ops;
		delete planeDataConfig;
	}


	~SimulationInstance()
	{
		if (simPlaneData != nullptr)
		{
			std::cout << "Error compared to initial condition" << std::endl;
			std::cout << "Lmax error: " << (simPlaneData->prog_h_t0-simPlaneData->prog_h).toPhys().physical_reduce_max_abs() << std::endl;
			std::cout << "RMS error: " << (simPlaneData->prog_h_t0-simPlaneData->prog_h).toPhys().physical_reduce_rms() << std::endl;
		}

		clear();
	}



#endif

#if 1


	void run_timestep()
	{
		if (timestepControl->current_simulation_time + timestepControl->current_timestep_size > timestepControl->max_simulation_time)
			timestepControl->current_timestep_size = timestepControl->max_simulation_time - timestepControl->current_simulation_time;

		timeSteppers.master->run_timestep(
				simPlaneData->prog_h, simPlaneData->prog_u, simPlaneData->prog_v,
				timestepControl->current_timestep_size,
				timestepControl->current_simulation_time
			);

		double dt = timestepControl->current_timestep_size;

		// advance in time
		timestepControl->current_simulation_time += dt;
		timestepControl->current_timestep_nr++;

		if (misc->verbosity > 2)
			std::cout << timestepControl->current_timestep_nr << ": " << timestepControl->current_simulation_time/(60*60*24.0) << std::endl;

		max_error_h0 = (simPlaneData->prog_h-simPlaneData->prog_h_t0).toPhys().physical_reduce_max_abs();
	}



	void compute_error()
	{
	}



	bool should_quit()
	{
		if (timestepControl->max_timesteps_nr != -1 && timestepControl->max_timesteps_nr <= timestepControl->current_timestep_nr)
			return true;

		double diff = std::abs(timestepControl->max_simulation_time - timestepControl->current_simulation_time);

		if (	timestepControl->max_simulation_time != -1 &&
				(
						timestepControl->max_simulation_time <= timestepControl->current_simulation_time	||
						diff/timestepControl->max_simulation_time < 1e-11	// avoid numerical issues in time stepping if current time step is 1e-14 smaller than max time step
				)
			)
			return true;

		return false;
	}


#if SWEET_GUI
	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(int i_num_iterations)
	{
		if (timestepControl->run_simulation_timesteps)
			for (int i = 0; i < i_num_iterations && !should_quit(); i++)
				run_timestep();

		compute_error();
	}


	void vis_get_vis_data_array(
			const PlaneData_Physical **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive_id,
			void **o_bogus_data,
			double *o_viz_min,
			double *o_viz_max,
			bool *viz_reset
	)
	{
		*o_render_primitive_id = render_primitive_id;

		int id = misc->vis_id % 3;
		switch (id)
		{
		case 0:
			sweet::convert(simPlaneData->prog_h, viz_plane_data);
			break;

		case 1:
			sweet::convert(simPlaneData->prog_u, viz_plane_data);
			break;

		case 2:
			sweet::convert(simPlaneData->prog_v, viz_plane_data);
			break;
		}

		*o_dataArray = &viz_plane_data;
		*o_aspect_ratio = 1;
	}



	const char* vis_get_status_string()
	{
		const char* description = "";
		int id = misc->vis_id % 3;

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
				"Time: %f (%.2f d), k: %i, dt: %.3e, Vis: %s, TMass: %.6e, TEnergy: %.6e, PotEnstrophy: %.6e, MaxVal: %.6e, MinVal: %.6e ",
#if SWEET_MPI
				-1,	// TODO: mpi_rank,
#endif
				timestepControl->current_simulation_time,
				timestepControl->current_simulation_time/(60.0*60.0*24.0),
				timestepControl->current_timestep_nr,
				timestepControl->current_timestep_size,
				description,
				diagnostics->total_mass,
				diagnostics->total_energy,
				diagnostics->total_potential_enstrophy,
				viz_plane_data.physical_reduce_max(),
				viz_plane_data.physical_reduce_min()
		);

		return title_string;
	}

	void vis_pause()
	{
		timestepControl->run_simulation_timesteps = !timestepControl->run_simulation_timesteps;
	}


	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			misc->vis_id++;
			break;

		case 'V':
			misc->vis_id--;
			break;

		case 'b':
			render_primitive_id = (render_primitive_id + 1) % 2;
			break;
		}
	}
#endif

#endif
};



int main(int i_argc, char *i_argv[])
{
	SimulationInstance simulation;
	if (simulation.error.exists())
	{
		std::cerr << "ERROR: " << simulation.error.get() << std::endl;
		return EXIT_FAILURE;
	}

	simulation.setup(i_argc, i_argv);
	if (simulation.error.exists())
	{
		std::cerr << "ERROR: " << simulation.error.get() << std::endl;
		return EXIT_FAILURE;
	}


#if SWEET_GUI
	if (simulation.gui_active)
	{
		std::cout << "GUI" << std::endl;
		//planeDataConfigInstance.setupAutoSpectralSpace(simVars.disc.space_res_physical, simVars.misc.reuse_spectral_transformation_plans);

		VisSweet<SimulationInstance> visSweet(&simulation);
		std::cout << "Max error h0: "<< simulation.max_error_h0 << std::endl;
	}
	else
#endif
	{
#if 0
		while (!simulation.should_quit())
		{
			simulation.run_timestep();
		}
#endif
	}

	std::cout << "FIN" << std::endl;
	return 0;
}
