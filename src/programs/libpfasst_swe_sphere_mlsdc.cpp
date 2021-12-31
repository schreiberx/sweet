/*
 * Author: Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst_swe_sphere_mlsdc
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */


#include <mpi.h>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>
#include "libpfasst_interface/LevelSingleton.hpp"
#include "libpfasst_swe_sphere_mlsdc/SphereDataCtx.hpp"
#include "swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <sweet/SWEETError.hpp>

#define WITH_MPI

SimulationVariables simVars;
std::vector<LevelSingleton> levelSingletons;

extern "C"
{
/* Driver function for pfasst control */
void fmain (SphereDataCtx* pd_ctx,
		const int*     nlevels,
		const int*     niters,
		const int*     nsweeps_coarse,
		const int      nnodes[],
		const char*    qtype_name,
		const int*     qtype_name_len,
		const int*     use_rk_stepper, // 1 means true, 0 means false
		const int*     nfields,
		const int      nvars_per_field[],
		double*        t_max,
		double*        dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
	MPI_Init(&i_argc, &i_argv);

	std::cout << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "WARNING: This program is NOT ready for use!" << std::endl;
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++" << std::endl;

	// input parameter names (specific ones for this program)
	const char *bogus_var_names[] = {
			"compute-error",
			nullptr
	};

	// set output time scale to hours
	simVars.iodata.output_time_scale = 1.0/(60.0*60.0);

	// default values for specific input (for general input see SimulationVariables.hpp)
	simVars.bogus.var[0] = 1;

	// Help menu
	if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
	{
		std::cout << "--compute-error [0/1]Output errors (if available, default: 1)" << std::endl;
		return -1;
	}

	// define the number of levels and SDC nodes for each level
	// note: level #nlevels-1 is the finest, level #0 is the coarsest

	int nnodes[simVars.libpfasst.nlevels];
	nnodes[simVars.libpfasst.nlevels-1] = simVars.libpfasst.nnodes; // finest level

	switch (simVars.libpfasst.nlevels)
	{
	// One level (nothing to do)
	case 1: {
		break;
	}
	// Two levels
	case 2: {
		if (simVars.libpfasst.nnodes == 3)
			nnodes[0] = 2;
		else if (simVars.libpfasst.nnodes == 5 ||
				simVars.libpfasst.nnodes == 9)
			nnodes[0] = 3;
		else if (simVars.libpfasst.nnodes == 7)
			nnodes[0] = 4; // for rk_stepper
		else
			SWEETError("With 2 levels, the number of SDC nodes on the fine level must be either 3, 5, or 9");
		break;
	}
	// Three levels
	case 3: {
		if (simVars.libpfasst.nnodes == 9)
		{
			nnodes[0] = 3;
			nnodes[1] = 5;
		}
		else if (simVars.libpfasst.nnodes == 5)
		{
			nnodes[0] = 2;
			nnodes[1] = 3;
		}
		else if (simVars.libpfasst.nnodes == 3)
		{
			nnodes[0] = 3;
			nnodes[1] = 3;
		}
		else
			SWEETError("With 3 levels, the number of SDC nodes on the fine level must be either 5, or 9");
		break;
	}
	// All other cases not supported yet
	default:
		SWEETError("Only 1, 2, or 3 levels are currently supported");
	}

	// setup the LevelSingletons for all levels
	// note: level #nlevels-1 is the finest, level #0 is the coarsest

	levelSingletons.resize(simVars.libpfasst.nlevels);

	// setup the finest level singleton
	const int fineLevelId = simVars.libpfasst.nlevels-1;


	levelSingletons[fineLevelId].level = fineLevelId;

	// setup data configuration in fine level

	levelSingletons[fineLevelId].dataConfig.setupAuto(
			simVars.disc.space_res_physical,
			simVars.disc.space_res_spectral,
			simVars.misc.reuse_spectral_transformation_plans,
			simVars.misc.verbosity
	);
	std::cout << "SPH config string: " << levelSingletons[fineLevelId].dataConfig.getConfigInformationString() << std::endl;

	int res_physical_nodealiasing[2] = {
			2*(simVars.disc.space_res_spectral[0]+1),
			simVars.disc.space_res_spectral[1]+2
	};

	levelSingletons[fineLevelId].dataConfigNoDealiasing.setupAuto(
			res_physical_nodealiasing,
			simVars.disc.space_res_spectral,
			simVars.misc.reuse_spectral_transformation_plans
	);

	// setup data operators in fine level

	levelSingletons[fineLevelId].op.setup(
			&(levelSingletons[fineLevelId].dataConfig),
			&(simVars.sim)
	);
	levelSingletons[fineLevelId].opNoDealiasing.setup(
			&(levelSingletons[fineLevelId].dataConfigNoDealiasing),
			&(simVars.sim)
	);

	// define the number of modes for the coarser levels
	for (int i = 1; i < simVars.libpfasst.nlevels; i++)
	{
		const int thisLevelId = simVars.libpfasst.nlevels-1-i;
		levelSingletons[thisLevelId].level = thisLevelId;

        // compute "additional" modes (negative because we're coarsening)
		auto additional_modes_lat = 1 - std::ceil(simVars.disc.space_res_spectral[0]*pow(simVars.libpfasst.coarsening_multiplier,i));
        auto additional_modes_lon = 1 - std::ceil(simVars.disc.space_res_spectral[1]*pow(simVars.libpfasst.coarsening_multiplier,i));
        // setup data configuration at this level
		levelSingletons[thisLevelId].dataConfig.setupAdditionalModes(
				&(levelSingletons[simVars.libpfasst.nlevels-i].dataConfig),
				additional_modes_lat,
				additional_modes_lon,
				simVars.misc.reuse_spectral_transformation_plans
		);

		// setup data operators at this level

		levelSingletons[thisLevelId].op.setup(
				&(levelSingletons[thisLevelId].dataConfig),
				&(simVars.sim)
		);
	}

	// define the SWEET parameters

	const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
	int nvars_per_field[simVars.libpfasst.nlevels];
	for (int i = 0; i < simVars.libpfasst.nlevels; ++i)
		nvars_per_field[i] = 2*levelSingletons[i].dataConfig.spectral_array_data_number_of_elements;  // number of degrees of freedom per vector field

	// initialize the topography before instantiating the SphereDataCtx object
	if (simVars.benchmark.benchmark_name == "flow_over_mountain")
	{

		// create h_topo with the configuration at the finest level
		simVars.benchmark.h_topo = SphereData_Physical(&(levelSingletons[simVars.libpfasst.nlevels-1].dataConfig));

		// initialize the topography
		(levelSingletons[simVars.libpfasst.nlevels-1].benchmarks).master->setup_topography();
	}

	// instantiate the SphereDataCtx object
	SphereDataCtx* pd_ctx = new SphereDataCtx(
			&simVars,
			&levelSingletons,
			nnodes
	);

	// get the C string length (needed by Fortran...)
	int string_length = simVars.libpfasst.nodes_type.size();

	// flag for the RK stepper
	const int rk_stepper_flag = (simVars.libpfasst.use_rk_stepper)
                            		? 1
                            				: 0;

	// call LibPFASST to advance in time
	fmain(
			pd_ctx,                                       // user defined context
			&simVars.libpfasst.nlevels,                   // number of SDC levels
			&simVars.libpfasst.niters,                    // number of SDC iterations
			&simVars.libpfasst.nsweeps_coarse,            // number of SDC sweeps on coarse level
			nnodes,                                       // number of SDC nodes
			(simVars.libpfasst.nodes_type).c_str(),       // type of nodes
			&string_length,                               // length of (simVars.libpfasst.nodes_type).c_str()
			&rk_stepper_flag,                             // flag for the RK stepper => 1 means true, 0 means false
			&nfields,                                     // number of vector fields
			nvars_per_field,                              // number of dofs per vector field
			&(simVars.timecontrol.max_simulation_time),   // simulation time
			&(simVars.timecontrol.current_timestep_size)  // time step size
	);

	// release the memory
	delete pd_ctx;

	MPI_Finalize();
}

