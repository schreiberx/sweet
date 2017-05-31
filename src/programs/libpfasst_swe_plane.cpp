/*
 * main.cpp
 *
 * PFASST SWE on the plane implementation
 *
 *  Created on: 30 Nov 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#include <sweet/FatalError.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/plane/PlaneDataTimesteppingRK.hpp>
#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/plane/PlaneDiagnostics.hpp>
#include <benchmarks_plane/PlaneBenchmarksCombined.hpp>

#include "libpfasst_swe_plane/LevelSingleton.hpp"
#include "libpfasst_swe_plane/PlaneDataCtx.hpp"
#include <mpi.h>

#define WITH_MPI

#if !SWEET_PFASST_CPP
#	error "Use --pfasst-cpp=enable for this file"
#endif

/**
 * Global variables which are shared between everything
 */

bool   param_use_coriolis_formulation       = true;
bool   param_compute_error                  = false;
double param_geostr_balance_freq_multiplier = 1.0;
int    param_timestepping_mode              = 0.0;
int    param_max_levels                     = 0;

SimulationVariables simVars;
std::vector<LevelSingleton> levelSingletons;

extern "C"
{
  /* Driver function for pfasst control */
  void fmain (PlaneDataCtx *pd_ctx, int* num_levs, std::size_t* nvars);
}

int main(int i_argc, char *i_argv[])
{
  MPI_Init(&i_argc, &i_argv);

  // input parameter names (specific ones for this program)
  const char *bogus_var_names[] = {
    "rexi-use-coriolis-formulation",
    "compute-error",
    nullptr
  };

  // default values for specific input (for general input see SimulationVariables.hpp)
  simVars.bogus.var[0] = 1;
  simVars.bogus.var[1] = 1;
  
  // Help menu
  if (!simVars.setupFromMainParameters(i_argc, i_argv, bogus_var_names))
    {
      std::cout << "	--compute-error [0/1]	Output errors (if available, default: 1)" << std::endl;
      std::cout << "	--rexi-use-coriolis-formulation [0/1]	Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
      return -1;
    }
  
  param_use_coriolis_formulation = simVars.bogus.var[0];
  assert (param_use_coriolis_formulation == 0 || param_use_coriolis_formulation == 1);
  param_compute_error = simVars.bogus.var[1];

  /**********************************************************
   * SETUP the LevelSingletons for all levels
   **********************************************************/
  
  // currently only one level
  param_max_levels = 1;
  
  // setup the level singleton
  levelSingletons.resize(param_max_levels);

  levelSingletons[0].dataConfig.setupAutoPhysicalSpace(
						       simVars.disc.res_spectral[0],
						       simVars.disc.res_spectral[1]
						       );

  levelSingletons[0].level = 0;
  
  levelSingletons[0].op.setup(
			      &(levelSingletons[0].dataConfig),
			      simVars.sim.domain_size,
			      simVars.disc.use_spectral_basis_diffs
			      );
  
  // this is not doing anything for now
  for (int i = 1; i < param_max_levels; i++)
    {
      levelSingletons[i].dataConfig.setupAdditionalModes(
							 &(levelSingletons[i-1].dataConfig),
							 -1,
							 -1
							 );
    }

  /**********************************************************
   * SETUP top level data
   **********************************************************/

  // instantiate the PlaneDataCtx object 
  PlaneDataCtx* pd_ctx = new PlaneDataCtx(
					  &simVars,
					  &levelSingletons
					  );
  
  // run a timestep with libpfasst and check that the solution matches the sweet integrators
  fmain(
	pd_ctx,
	&param_max_levels,
	&(levelSingletons[0].dataConfig.physical_array_data_number_of_elements)
	); 
  
  // release the memory
  delete pd_ctx;

  MPI_Finalize();
}

