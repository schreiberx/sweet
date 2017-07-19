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

SimulationVariables simVars;
std::vector<LevelSingleton> levelSingletons;

extern "C"
{
  /* Driver function for pfasst control */
  void fmain (PlaneDataCtx *pd_ctx, 
	      const int* nlevs, const int* niters, const int* nnodes, 
	      const int* nfields, const size_t* nvars_per_field, 
	      double* t_max, double* dt);
}

/**
 * Main function launching LibPFASST
 */

int main(int i_argc, char *i_argv[])
{
  MPI_Init(&i_argc, &i_argv);

  simVars.timecontrol.current_timestep_size = - simVars.sim.CFL; 
  simVars.outputConfig();

  // define the LibPFASST parameters (later implemented as command line args)

  const int nlevels  = 2;                       // number of SDC levels
  const int niters   = 2;                       // number of SDC iterations
  const int nnodes[] = {3, 5};                  // number of SDC nodes
  const double coarsening_multiplier[] = { 1 }; // spatial coarsening ratio for the levels



  // setup the LevelSingletons for all levels
  // note: level #nlevels-1 is the finest, level #0 is the coarsest

  // setup the finest level singleton
  levelSingletons.resize(nlevels);

  levelSingletons[nlevels-1].dataConfig.setupAutoPhysicalSpace(
							       simVars.disc.res_spectral[0],
							       simVars.disc.res_spectral[1]
							       );

  levelSingletons[nlevels-1].level = nlevels-1;
  
  levelSingletons[nlevels-1].op.setup(
				      &(levelSingletons[nlevels-1].dataConfig),
				      simVars.sim.domain_size,
				      simVars.disc.use_spectral_basis_diffs
				      );
  
  // define the number of modes for the coarser levels
  for (int i = 1; i < nlevels; i++)
    {
      levelSingletons[nlevels-1-i].dataConfig.setupAdditionalModes(
								   &(levelSingletons[nlevels-i].dataConfig),
								   -simVars.disc.res_spectral[0]*coarsening_multiplier[i-1],
								   -simVars.disc.res_spectral[1]*coarsening_multiplier[i-1]
								   );
      
      levelSingletons[nlevels-1-i].level = nlevels-1-i;
  
      levelSingletons[nlevels-1-i].op.setup(
					    &(levelSingletons[nlevels-1-i].dataConfig),
					    simVars.sim.domain_size,
					    simVars.disc.use_spectral_basis_diffs
					    );
    }



  // define the SWEET parameters

  const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
  const size_t nvars_per_field[] = { levelSingletons[0].dataConfig.physical_array_data_number_of_elements, // number of degrees of freedom per vector field
				     levelSingletons[1].dataConfig.physical_array_data_number_of_elements };



  // instantiate the PlaneDataCtx object 

  PlaneDataCtx* pd_ctx = new PlaneDataCtx(
  					  &simVars,
  					  &levelSingletons
  					  );
  
  // output the info for the levels
  for (int i = 0; i < nlevels; i++)
    levelSingletons[nlevels-1-i].dataConfig.printInformation();
  


  // call LibPFASST to advance in time
  fmain(
	pd_ctx,                                       // user defined context
	&nlevels,                                     // number of SDC levels
	&niters,                                      // number of SDC iterations
	&nnodes[0],                                   // number of SDC nodes 
	&nfields,                                     // number of vector fields
	&nvars_per_field[0],                          // number of dofs per vector field
	&(simVars.timecontrol.max_simulation_time),   // simulation time
	&(simVars.timecontrol.current_timestep_size)  // time step size
  	); 

  // release the memory
  delete pd_ctx;

  MPI_Finalize();
}

