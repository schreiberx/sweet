/*
 * Author: Valentina Schüller & Francois Hamon & Martin Schreiber <SchreiberX@gmail.com>
 * MULE_COMPILE_FILES_AND_DIRS: src/programs/libpfasst/pde_sweSphere_imex_sdc/
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

#include <mpi.h>
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>
#include <sweet/core/shacks/ShackProgArgDictionary.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>
#include "../pde_sweSphere/ShackPDESWESphere.hpp"
#include "../pde_sweSphere/benchmarks/ShackPDESWESphereBenchmarks.hpp"
#include "ShackLibPFASST.hpp"

#include <sweet/core/SWEETError.hpp>
#include "interface/LevelSingleton.hpp"
#include "pde_sweSphere_imex_sdc/SphereDataCtxSDC.hpp"
#include "../pde_sweSphere/PDESWESphere_Diagnostics.hpp"

#define WITH_MPI

extern "C"
{
    /* Driver function for pfasst control */
    void fmain(SphereDataCtxSDC *pd_ctx,
               const int *niters,
               const int *nsweeps,
               const int *nnodes,
               const char *qtype_name,
               const int *qtype_name_len,
               const int *use_rk_stepper, // 1 means true, 0 means false
               const int *nfields,
               const int *nvars_per_field,
               double *t_max,
               double *dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
    MPI_Init(&i_argc, &i_argv);

    sweet::ShackProgArgDictionary shackProgArgDict(i_argc, i_argv);
    shackProgArgDict.setup();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    sweet::ShackSphereDataOps *shackSphereDataOps = shackProgArgDict.getAutoRegistration<sweet::ShackSphereDataOps>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    ShackPDESWESphere *shackPDESWESphere = shackProgArgDict.getAutoRegistration<ShackPDESWESphere>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    ShackLibPFASST *shackLibPFASST = shackProgArgDict.getAutoRegistration<ShackLibPFASST>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    // sweet::ShackTimestepControl *shackTimestepControl =
    shackProgArgDict.getAutoRegistration<sweet::ShackTimestepControl>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    // ShackPDESWESphereBenchmarks *shackBenchmarks =
    shackProgArgDict.getAutoRegistration<ShackPDESWESphereBenchmarks>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    // sweet::ShackIOData *shackIOData =
    shackProgArgDict.getAutoRegistration<sweet::ShackIOData>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    // ShackPDESWESphereTimeDiscretization *shackTimeDisc =
    shackProgArgDict.getAutoRegistration<ShackPDESWESphereTimeDiscretization>();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    shackProgArgDict.processProgramArguments();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    shackProgArgDict.checkAllArgumentsProcessed();
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(shackProgArgDict);

    shackProgArgDict.printShackData();

    sweet::SphereData_Config sphereData_Config;
    sphereData_Config.setupAuto(shackSphereDataOps);
    ERROR_CHECK_WITH_PRINT_AND_COND_RETURN_EXITCODE(sphereData_Config);

    LevelSingleton levelSingleton;

    if (shackLibPFASST->nlevels != 1)
    {
        SWEETError("For SDC, nlevels has to be equal to 1");
    }

    shackLibPFASST->postprocess_nsweeps();

    // set up levelSingleton
    levelSingleton.level = 0;
    levelSingleton.sphereDataConfig.setupAuto(shackSphereDataOps);
    std::cout << "SPH config string: " << levelSingleton.sphereDataConfig.getConfigInformationString() << std::endl;

    // setup data operators
    levelSingleton.ops.setup(
        &(levelSingleton.sphereDataConfig),
        shackSphereDataOps);

    // define the SWEET parameters

    const int nfields = 3; // number of vector fields (here, height and two horizontal velocities)
    int nvars_per_field;
    nvars_per_field = 2 * levelSingleton.sphereDataConfig.spectral_array_data_number_of_elements; // number of degrees of freedom per vector field

    // initialize the topography before instantiating the SphereDataCtxSDC object
    /*if (shackDict.benchmark.benchmark_name == "flow_over_mountain")
    {
        // create h_topo with the configuration at the finest level
        shackDict.benchmark.h_topo = SphereData_Physical(&(levelSingleton.sphereDataConfig));

        // initialize the topography
        levelSingleton.benchmarks.master->setup_topography();
    }*/

    // instantiate the SphereDataCtxSDC object
    int nnodes[1];
    nnodes[0] = shackLibPFASST->nnodes;
    SphereDataCtxSDC pd_ctx;

    pd_ctx.shackRegistration(&shackProgArgDict);

    pd_ctx.setup(
        &shackProgArgDict,
        &levelSingleton,
        nnodes);

    // get the C string length (needed by Fortran...)
    int string_length = shackLibPFASST->nodes_type.size();

    // flag for the RK stepper
    const int rk_stepper_flag = (shackLibPFASST->use_rk_stepper) ? 1 : 0;

    // call LibPFASST to advance in time
    fmain(
        &pd_ctx,                                             // user defined context
        &shackLibPFASST->niters,                             // number of SDC iterations
        &shackLibPFASST->nsweeps[0],                         // number of SDC sweeps on coarse level
        &nnodes[0],                                          // number of SDC nodes
        (shackLibPFASST->nodes_type).c_str(),                // type of nodes
        &string_length,                                      // length of (shackLibPFASST->nodes_type).c_str()
        &rk_stepper_flag,                                    // flag for the RK stepper => 1 means true, 0 means false
        &nfields,                                            // number of vector fields
        &nvars_per_field,                                    // number of dofs per vector field
        &(pd_ctx.shackTimestepControl->max_simulation_time), // simulation time
        &(pd_ctx.shackTimestepControl->current_timestepSize));

    MPI_Finalize();
}
