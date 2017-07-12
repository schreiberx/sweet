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

// #if !SWEET_PFASST_CPP
// #	error "Use --pfasst-cpp=enable for this file"
// #endif

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
  void fmain (PlaneDataCtx *pd_ctx, int* num_levs, double* t_max, double* dt, std::size_t* nvars);
}

/**
 * Main function launching LibPFASST
 */
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

  simVars.timecontrol.current_timestep_size = - simVars.sim.CFL; 
  simVars.outputConfig();

  /**********************************************************
   * SETUP the LevelSingletons for all levels
   **********************************************************/
  
  // note: level #param_max_levels-1 is the finest, level #0 is the coarsest

  // define the number of levels
  param_max_levels = 1;
  
  // setup the level singleton
  levelSingletons.resize(param_max_levels);

  levelSingletons[param_max_levels-1].dataConfig.setupAutoPhysicalSpace(
									simVars.disc.res_spectral[0],
									simVars.disc.res_spectral[1]
									);

  levelSingletons[param_max_levels-1].level = param_max_levels-1;
  
  levelSingletons[param_max_levels-1].op.setup(
					       &(levelSingletons[param_max_levels-1].dataConfig),
					       simVars.sim.domain_size,
					       simVars.disc.use_spectral_basis_diffs
					       );
  
  // define the number of modes for the coarser levels
  int denominator = 2;
  for (int i = 1; i < param_max_levels; i++)
    {
      levelSingletons[param_max_levels-1-i].dataConfig.setupAdditionalModes(
									    &(levelSingletons[param_max_levels-i].dataConfig),
									    -simVars.disc.res_spectral[0]/denominator,
									    -simVars.disc.res_spectral[1]/denominator
									    );

      levelSingletons[param_max_levels-1-i].level = param_max_levels-1-i;
  
      levelSingletons[param_max_levels-1-i].op.setup(
						     &(levelSingletons[param_max_levels-1-i].dataConfig),
						     simVars.sim.domain_size,
						     simVars.disc.use_spectral_basis_diffs
						     );

      denominator *= 2;

    }

  /**********************************************************
   * SETUP top level data
   **********************************************************/

  // instantiate the PlaneDataCtx object 
  PlaneDataCtx* pd_ctx = new PlaneDataCtx(
  					  &simVars,
  					  &levelSingletons
  					  );
  
  // output the info for the fine level
  levelSingletons[param_max_levels-1].dataConfig.printInformation();
  
  // run a timestep with libpfasst and check that the solution matches the sweet integrators
  fmain(
	pd_ctx,
	&param_max_levels,
	&(simVars.timecontrol.max_simulation_time),
	&(simVars.timecontrol.current_timestep_size),
	&(levelSingletons[param_max_levels-1].dataConfig.physical_array_data_number_of_elements)
  	); 

  // release the memory
  delete pd_ctx;

  MPI_Finalize();
}

