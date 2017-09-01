
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
	      const int* nlevels, const int* niters, const int nnodes[], const char* qtype_name, const int* qtype_name_len,
	      const int* nfields, const int nvars_per_field[], 
	      double* t_max, double* dt);
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
      std::cout << "--compute-error [0/1]Output errors (if available, default: 1)" << std::endl;
      std::cout << "--rexi-use-coriolis-formulation [0/1]Use Coriolisincluding  solver for REXI (default: 1)" << std::endl;
      return -1;
    }

  simVars.timecontrol.current_timestep_size = - simVars.sim.CFL; 
  simVars.outputConfig();



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
	if (simVars.libpfasst.nnodes == 3 ||
	    simVars.libpfasst.nnodes == 5 || 
	    simVars.libpfasst.nnodes == 9)
	  nnodes[0] = 3; 
	else 
	  FatalError("With 2 levels, the number of SDC nodes on the fine level must be either 3, 5, or 9");
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
	else 
	  FatalError("With 3 levels, the number of SDC nodes on the fine level must be either 5, or 9");
	break;
      }
      // All other cases not supported yet
      default:
	FatalError("Only 1, 2, or 3 levels are currently supported");
    }
  
  // setup the LevelSingletons for all levels
  // note: level #nlevels-1 is the finest, level #0 is the coarsest



  // setup the finest level singleton
  levelSingletons.resize(simVars.libpfasst.nlevels);

  levelSingletons[simVars.libpfasst.nlevels-1].dataConfig.setupAutoPhysicalSpace(
										 simVars.disc.res_spectral[0],
										 simVars.disc.res_spectral[1],
										 &(simVars.disc.res_physical[0]),
										 &(simVars.disc.res_physical[1])
										 );

  levelSingletons[simVars.libpfasst.nlevels-1].level = simVars.libpfasst.nlevels-1;
  
  levelSingletons[simVars.libpfasst.nlevels-1].op.setup(
							&(levelSingletons[simVars.libpfasst.nlevels-1].dataConfig),
							simVars.sim.domain_size,
							simVars.disc.use_spectral_basis_diffs
							);
  
  // define the number of modes for the coarser levels
  for (int i = 1; i < simVars.libpfasst.nlevels; i++)
    {
      levelSingletons[simVars.libpfasst.nlevels-1-i].dataConfig.setupAdditionalModes(
										     &(levelSingletons[simVars.libpfasst.nlevels-i].dataConfig),
										     -simVars.disc.res_spectral[0]*pow(simVars.libpfasst.coarsening_multiplier,i),
										     -simVars.disc.res_spectral[1]*pow(simVars.libpfasst.coarsening_multiplier,i)
										     );
      
      levelSingletons[simVars.libpfasst.nlevels-1-i].level = simVars.libpfasst.nlevels-1-i;
  
      levelSingletons[simVars.libpfasst.nlevels-1-i].op.setup(
							      &(levelSingletons[simVars.libpfasst.nlevels-1-i].dataConfig),
							      simVars.sim.domain_size,
							      simVars.disc.use_spectral_basis_diffs
							      );
    }



  // define the SWEET parameters

  const int nfields = 3;  // number of vector fields (here, height and two horizontal velocities)
  int nvars_per_field[simVars.libpfasst.nlevels];
  for (int i = 0; i < simVars.libpfasst.nlevels; ++i) 
    nvars_per_field[i] = levelSingletons[i].dataConfig.physical_array_data_number_of_elements;  // number of degrees of freedom per vector field

  // instantiate the PlaneDataCtx object 
  PlaneDataCtx* pd_ctx = new PlaneDataCtx(
  					  &simVars,
  					  &levelSingletons,
					  nnodes
  					  );
  // output the info for the levels
  for (int i = 0; i < simVars.libpfasst.nlevels; i++)
    levelSingletons[simVars.libpfasst.nlevels-1-i].dataConfig.printInformation();
  
  // get the C string length (needed by Fortran...)
  int string_length = simVars.libpfasst.nodes_type.size();

  // call LibPFASST to advance in time
  fmain(
	pd_ctx,                                       // user defined context
	&simVars.libpfasst.nlevels,                   // number of SDC levels
	&simVars.libpfasst.niters,                    // number of SDC iterations
	nnodes,                                       // number of SDC nodes 
	(simVars.libpfasst.nodes_type).c_str(),       // type of nodes
	&string_length,                               // length of (simVars.libpfasst.nodes_type).c_str()
	&nfields,                                     // number of vector fields
	nvars_per_field,                              // number of dofs per vector field
	&(simVars.timecontrol.max_simulation_time),   // simulation time
	&(simVars.timecontrol.current_timestep_size)  // time step size
  	); 

  // release the memory
  delete pd_ctx;

  MPI_Finalize();
}

