/*
 * Author: Valentina Schüller & Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst_swe_sphere_expl_sdc
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/swe_sphere_timeintegrators/
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 */


#include <mpi.h>
#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <sweet/SimulationVariables.hpp>
#include "libpfasst_interface/LevelSingleton.hpp"
#include "libpfasst_swe_sphere_expl_sdc/SphereDataCtxSDC.hpp"
#include "swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <sweet/SWEETError.hpp>

#define WITH_MPI

SimulationVariables simVars;
LevelSingleton levelSingleton;

extern "C"
{
/* Driver function for pfasst control */
void fmain (SphereDataCtxSDC* pd_ctx,
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

	if (simVars.libpfasst.nlevels != 1)
	{
		SWEETError("For SDC, nlevels has to be equal to 1");
	}

	// set up level of levelSingleton
	levelSingleton.level = 0;

	// setup data configuration in fine level

	levelSingleton.dataConfig.setupAuto(
			simVars.disc.space_res_physical,
			simVars.disc.space_res_spectral,
			simVars.misc.reuse_spectral_transformation_plans,
			simVars.misc.verbosity
	);
	std::cout << "SPH config string: " << levelSingleton.dataConfig.getConfigInformationString() << std::endl;

	int res_physical_nodealiasing[2] = {
			2*(simVars.disc.space_res_spectral[0]+1),
			simVars.disc.space_res_spectral[1]+2
	};

	levelSingleton.dataConfigNoDealiasing.setupAuto(
			res_physical_nodealiasing,
			simVars.disc.space_res_spectral,
			simVars.misc.reuse_spectral_transformation_plans
	);

	// setup data operators
	levelSingleton.op.setup(
			&(levelSingleton.dataConfig),
			&(simVars.sim)
	);
	levelSingleton.opNoDealiasing.setup(
			&(levelSingleton.dataConfigNoDealiasing),
			&(simVars.sim)
	);


	// define the SWEET parameters

	const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
	int nvars_per_field[1];
	nvars_per_field[0] = 2 * levelSingleton.dataConfig.spectral_array_data_number_of_elements;  // number of degrees of freedom per vector field

	// initialize the topography before instantiating the SphereDataCtxSDC object
	if (simVars.benchmark.benchmark_name == "flow_over_mountain")
	{
		// create h_topo with the configuration at the finest level
		simVars.benchmark.h_topo = SphereData_Physical(&(levelSingleton.dataConfig));

		// initialize the topography
		levelSingleton.benchmarks.master->setup_topography();
	}

	// instantiate the SphereDataCtxSDC object
	SphereDataCtxSDC* pd_ctx = new SphereDataCtxSDC(
			&simVars,
			&levelSingleton,
			nnodes
	);

	// get the C string length (needed by Fortran...)
	int string_length = simVars.libpfasst.nodes_type.size();

	// flag for the RK stepper
	const int rk_stepper_flag = (simVars.libpfasst.use_rk_stepper) ? 1 : 0;

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

