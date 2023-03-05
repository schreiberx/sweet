/*
 * Author: Valentina Schüller & Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst_swe_sphere_imex_sdc
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 * MULE_SCONS_OPTIONS: --lapack=enable
 * MULE_SCONS_OPTIONS: --sweet-mpi=enable
 * MULE_SCONS_OPTIONS: --libpfasst=enable
 */


#include <mpi.h>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>
#include "libpfasst_interface/LevelSingleton.hpp"
#include "libpfasst_swe_sphere_imex_sdc/SphereDataCtxSDC.hpp"
#include "swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <sweet/core/SWEETError.hpp>

#define WITH_MPI

extern "C"
{
/* Driver function for pfasst control */
void fmain (SphereDataCtxSDC* pd_ctx,
		const int*     niters,
		const int*     nsweeps,
		const int*     nnodes,
		const char*    qtype_name,
		const int*     qtype_name_len,
		const int*     use_rk_stepper, // 1 means true, 0 means false
		const int*     nfields,
		const int*     nvars_per_field,
		double*        t_max,
		double*        dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
	MPI_Init(&i_argc, &i_argv);

	sweet::ShackDictionary shackDict;
	LevelSingleton levelSingleton;

	// input parameter names (specific ones for this program)
	const char *user_defined_prog_params[] = {
			"compute-errors",
			nullptr
	};

	// set output time scale to hours
	shackDict.iodata.output_time_scale = 1.0/(60.0*60.0);

	// default values for specific input (for general input see shacks/ShackDictionary.hpp)
	shackDict.user_defined.var[0] = 1;

	// Help menu
	if (!shackDict.setupFromMainParameters(i_argc, i_argv, user_defined_prog_params))
	{
		std::cout << "--compute-errors [0/1]Output errors (if available, default: 1)" << std::endl;
		return -1;
	}

	if (shackDict.libpfasst.nlevels != 1)
	{
		SWEETError("For SDC, nlevels has to be equal to 1");
	}

	shackDict.libpfasst.postprocess_nsweeps();

	// set up levelSingleton
	levelSingleton.level = 0;
	levelSingleton.dataConfig.setupAuto(
			shackDict.disc.space_res_physical,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackDict.misc.verbosity,
			shackDict.parallelization.num_threads_space
	);
	std::cout << "SPH config string: " << levelSingleton.dataConfig.getConfigInformationString() << std::endl;

	int res_physical_nodealiasing[2] = {
			2*(shackDict.disc.space_res_spectral[0]+1),
			shackDict.disc.space_res_spectral[1]+2
	};

	levelSingleton.dataConfigNoDealiasing.setupAuto(
			res_physical_nodealiasing,
			shackDict.disc.space_res_spectral,
			shackDict.misc.reuse_spectral_transformation_plans,
			shackDict.misc.verbosity,
			shackDict.parallelization.num_threads_space
	);

	// setup data operators
	levelSingleton.op.setup(
			&(levelSingleton.dataConfig),
			&(shackDict.sim)
	);
	levelSingleton.opNoDealiasing.setup(
			&(levelSingleton.dataConfigNoDealiasing),
			&(shackDict.sim)
	);


	// define the SWEET parameters

	const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
	int nvars_per_field;
	nvars_per_field = 2 * levelSingleton.dataConfig.spectral_array_data_number_of_elements;  // number of degrees of freedom per vector field

	// initialize the topography before instantiating the SphereDataCtxSDC object
	if (shackDict.benchmark.benchmark_name == "flow_over_mountain")
	{
		// create h_topo with the configuration at the finest level
		shackDict.benchmark.h_topo = SphereData_Physical(&(levelSingleton.dataConfig));

		// initialize the topography
		levelSingleton.benchmarks.master->setup_topography();
	}

	// instantiate the SphereDataCtxSDC object
	int nnodes[1];
	nnodes[0] = shackDict.libpfasst.nnodes;
	SphereDataCtxSDC* pd_ctx = new SphereDataCtxSDC(
			&shackDict,
			&levelSingleton,
			nnodes
	);

	// get the C string length (needed by Fortran...)
	int string_length = shackDict.libpfasst.nodes_type.size();

	// flag for the RK stepper
	const int rk_stepper_flag = (shackDict.libpfasst.use_rk_stepper) ? 1 : 0;

	// call LibPFASST to advance in time
	fmain(
			pd_ctx,                                       // user defined context
			&shackDict.libpfasst.niters,                    // number of SDC iterations
			&shackDict.libpfasst.nsweeps[0],            // number of SDC sweeps on coarse level
			&nnodes[0],                                       // number of SDC nodes
			(shackDict.libpfasst.nodes_type).c_str(),       // type of nodes
			&string_length,                               // length of (shackDict.libpfasst.nodes_type).c_str()
			&rk_stepper_flag,                             // flag for the RK stepper => 1 means true, 0 means false
			&nfields,                                     // number of vector fields
			&nvars_per_field,                              // number of dofs per vector field
			&(shackDict.timecontrol.max_simulation_time),   // simulation time
			&(shackDict.timecontrol.current_timestep_size)  // time step size
	);

	// release the memory
	delete pd_ctx;

	MPI_Finalize();
}

