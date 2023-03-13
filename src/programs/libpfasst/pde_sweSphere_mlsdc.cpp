/*
 * Author: Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 *
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst/pde_sweSphere_mlsdc/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/time/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/benchmarks/
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/pde_sweSphere/
 *
 * MULE_SCONS_OPTIONS: --sphere-spectral-space=enable
 * MULE_SCONS_OPTIONS: --fortran-source=enable
 * MULE_SCONS_OPTIONS: --lapack=enable
 * MULE_SCONS_OPTIONS: --sweet-mpi=enable
 * MULE_SCONS_OPTIONS: --libpfasst=enable
 */


#include <iostream>
#include <mpi.h>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "../pde_sweSphere/ShackPDESWESphere.hpp"
#include "../pde_sweSphere/benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "ShackLibPFASST.hpp"

#include "interface/LevelSingleton.hpp"
#include "pde_sweSphere_mlsdc/SphereDataCtx.hpp"
#include "../pde_sweSphere/PDESWESphere_BenchmarksCombined.hpp"
#include "../pde_sweSphere/PDESWESphere_Diagnostics.hpp"
#include "../pde_sweSphere/time/ShackPDESWESphereTimeDiscretization.hpp"
#include <sweet/core/SWEETError.hpp>



#define WITH_MPI

extern "C"
{
/* Driver function for pfasst control */
void fmain (
		SphereDataCtx* pd_ctx,
		const int*     nlevels,
		const int*     niters,
		const int      nsweeps[],
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


	sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
	shackProgArgDict.setup();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	ShackPDESWESphere *shackPDESWESphere = shackProgArgDict.getAutoRegistration<ShackPDESWESphere>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	ShackLibPFASST *shackLibPFASST = shackProgArgDict.getAutoRegistration<ShackLibPFASST>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	//sweet::ShackTimestepControl *shackTimestepControl =
			shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	//ShackPDESWESphereBenchmarks *shackBenchmarks =
			shackProgArgDict.getAutoRegistration<ShackPDESWESphereBenchmarks>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	//sweet::ShackIOData *shackIOData =
			shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	//ShackPDESWESphereTimeDiscretization *shackTimeDisc =
			shackProgArgDict.getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);


	shackProgArgDict.processProgramArguments();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.checkAllArgumentsProcessed();
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(shackProgArgDict);

	shackProgArgDict.printShackData();

	sweet::SphereData_Config sphereData_Config;
	sphereData_Config.setupAuto(shackSphereDataOps);
	ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXIT(sphereData_Config);


	std::vector<LevelSingleton> levelSingletons;

	// set output time scale to hours
	//shackIOData->output_time_scale = 1.0/(60.0*60.0);
	//shackIOData->output_time_scale_inv = 60.0*60.0;

	if ((shackPDESWESphere->viscosity > 0) || (shackPDESWESphere->viscosity_order != 2))
	{
		SWEETError("To apply viscosity, use the --libpfasst-u2/4/6/8 flags, not -u or -U!");
	}

	shackLibPFASST->postprocess_hyperviscosity();
	shackLibPFASST->postprocess_nsweeps();

	// define the number of levels and SDC nodes for each level
	// note: level #nlevels-1 is the finest, level #0 is the coarsest

	std::vector<int> nnodes;
	nnodes.resize(shackLibPFASST->nlevels);
	nnodes[shackLibPFASST->nlevels-1] = shackLibPFASST->nnodes; // finest level

	switch (shackLibPFASST->nlevels)
	{
	// One level (nothing to do)
	case 1: {
		break;
	}
	// Two levels
	case 2: {
		if (shackLibPFASST->nnodes == 3)
			nnodes[0] = 2;
		else if (shackLibPFASST->nnodes == 5 ||
				shackLibPFASST->nnodes == 9)
			nnodes[0] = 3;
		else if (shackLibPFASST->nnodes == 7)
			nnodes[0] = 4; // for rk_stepper
		else
			SWEETError("With 2 levels, the number of SDC nodes on the fine level must be either 3, 5, or 9");
		break;
	}
	// Three levels
	case 3: {
		if (shackLibPFASST->nnodes == 9)
		{
			nnodes[0] = 3;
			nnodes[1] = 5;
		}
		else if (shackLibPFASST->nnodes == 5)
		{
			nnodes[0] = 2;
			nnodes[1] = 3;
		}
		else if (shackLibPFASST->nnodes == 3)
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

	levelSingletons.resize(shackLibPFASST->nlevels);

	// setup the finest level singleton
	const int fineLevelId = shackLibPFASST->nlevels-1;


	levelSingletons[fineLevelId].level = fineLevelId;

	// setup data configuration in fine level

	levelSingletons[fineLevelId].sphereDataConfig.setupAuto(shackSphereDataOps);

	std::cout << "SPH config string: " << levelSingletons[fineLevelId].sphereDataConfig.getConfigInformationString() << std::endl;

	// setup data operators in fine level

	levelSingletons[fineLevelId].ops.setup(
			&(levelSingletons[fineLevelId].sphereDataConfig),
			shackSphereDataOps
	);

	// define the number of modes for the coarser levels
	for (int i = 1; i < shackLibPFASST->nlevels; i++)
	{
		const int thisLevelId = shackLibPFASST->nlevels-1-i;
		levelSingletons[thisLevelId].level = thisLevelId;

        // compute "additional" modes (negative because we're coarsening)
		// use 1 - alpha to compute what to take away (we want to have alpha^i * res modes on level n-1-i)
		double coarsener = 1 - shackLibPFASST->coarsening_multiplier;
		int additional_modes_lat = 1 - std::ceil(shackSphereDataOps->space_res_spectral[0]*pow(coarsener,i));
        int additional_modes_lon = 1 - std::ceil(shackSphereDataOps->space_res_spectral[1]*pow(coarsener,i));
        // setup data configuration at this level
		levelSingletons[thisLevelId].sphereDataConfig.setupAdditionalModes(
				&(levelSingletons[thisLevelId + 1].sphereDataConfig),
				additional_modes_lat,
				additional_modes_lon,
				shackSphereDataOps
		);

		// setup data operators at this level
		levelSingletons[thisLevelId].ops.setup(
				&(levelSingletons[thisLevelId].sphereDataConfig),
				shackSphereDataOps
		);
	}

	// define the SWEET parameters

	const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
	std::vector<int> nvars_per_field;
	nvars_per_field.resize(shackLibPFASST->nlevels);
	for (int i = 0; i < shackLibPFASST->nlevels; ++i)
	{
		// number of degrees of freedom per vector field
		nvars_per_field[i] = 2*levelSingletons[i].sphereDataConfig.spectral_array_data_number_of_elements; 
	}

	{
		// instantiate the SphereDataCtx object
		SphereDataCtx pd_ctx;

		pd_ctx.shackRegistration(&shackProgArgDict);

		pd_ctx.setup(
				&shackProgArgDict,
				&levelSingletons,
				&nnodes
		);

		// get the C string length (needed by Fortran...)
		int string_length = shackLibPFASST->nodes_type.size();

		// flag for the RK stepper
		const int rk_stepper_flag = shackLibPFASST->use_rk_stepper ? 1 : 0;

		// call LibPFASST to advance in time
		fmain(
				&pd_ctx,								// user defined context
				&shackLibPFASST->nlevels,			// number of SDC levels
				&shackLibPFASST->niters,			// number of SDC iterations
				shackLibPFASST->nsweeps.data(),		// number of SDC sweeps on coarse level
				nnodes.data(),						// number of SDC nodes
				(shackLibPFASST->nodes_type).c_str(),	// type of nodes
				&string_length,						// length of (shackLibPFASST->nodes_type).c_str()
				&rk_stepper_flag,					// flag for the RK stepper => 1 means true, 0 means false
				&nfields,							// number of vector fields
				nvars_per_field.data(),				// number of dofs per vector field
				&(pd_ctx.shackTimestepControl->max_simulation_time),   // simulation time
				&(pd_ctx.shackTimestepControl->current_timestep_size)  // time step size
		);
	}


	MPI_Finalize();
}

