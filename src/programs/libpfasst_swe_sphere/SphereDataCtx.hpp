#ifndef _SPHERE_DATA_CTX_HPP_
#define _SPHERE_DATA_CTX_HPP_

#include <sweet/sphere/SphereData_Spectral.hpp>
#include <sweet/sphere/SphereHelpers_Diagnostics.hpp>
#include <sweet/sphere/SphereOperators_SphereData.hpp>
#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"

#include "../swe_sphere_timeintegrators/SWE_Sphere_TS_ln_erk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the level singletons and the SimulationVariables object

class SphereDataCtx {

public:
  
  // Contructor
  SphereDataCtx(
               SimulationVariables *i_simVars,
	       std::vector<LevelSingleton> *i_singletons,
	       int* i_nnodes
	       ) 
    : simVars(i_simVars),
      levelSingletons(i_singletons)
  {
    int rank   = 0;
    int nprocs = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (!simVars) 
      SWEETError("SphereDataCtx: simVars pointer is NULL!");

    if (!levelSingletons) 
      SWEETError("SphereDataCtx: levelSingletons pointer is NULL!");

    // initialize the time steppers from SWEET

    timestepper_ln_erk = new SWE_Sphere_TS_ln_erk(
						  *simVars,
						  ((*levelSingletons)[levelSingletons->size()-1].op)
						  );
    
    int timestepping_order = i_simVars->disc.timestepping_order;
    if (timestepping_order == -1)
    {
      std::cout << "Timestepping Order not given! Choosing default = 4" << std::endl;
      timestepping_order = 4;
    }
    
    timestepper_ln_erk->setup(timestepping_order);
    
    // initialize the residuals
    residuals.resize(nprocs,std::vector<double>(0,0.));

    // initialize the diagnostics object
    sphereDiagnostics = new SphereHelpers_Diagnostics(
					      &((*levelSingletons)[levelSingletons->size()-1].dataConfig),
					      *simVars,
					      0
					      );
	
  }

  // Destructor
  ~SphereDataCtx() 
  {
    delete timestepper_ln_erk;
    delete sphereDiagnostics;
  }

  // Getter for the sphere data configuration at level i_level
  SphereData_Config* get_sphere_data_config(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].dataConfig);
  }

  // Getter for the sphere data configuration at level i_level
  BenchmarksSphereSWE* get_swe_benchmark(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].benchmarks);
  }

  // Getter for the sphere data configuration with no dealiasing at the fine level
  SphereData_Config* get_sphere_data_config_nodealiasing() const 
  {
    return &((*levelSingletons)[levelSingletons->size()-1].dataConfigNoDealiasing);
  }

  // Getter for the sphere data operators at level i_level
  SphereOperators_SphereData* get_sphere_operators(
					int i_level
					) const
  {
    return &((*levelSingletons)[i_level].op);
  }

  // Getter for the sphere data operators with no dealiasing at the fine level
  SphereOperators_SphereData* get_sphere_operators_nodealiasing() const
  {
    return &((*levelSingletons)[levelSingletons->size()-1].opNoDealiasing);
  }


  // Getter for the sphere diagnostics at the fine level
  SphereHelpers_Diagnostics* get_sphere_diagnostics() 
  {
    return sphereDiagnostics;
  }



  // Getter for the explicit timestepper
  SWE_Sphere_TS_ln_erk* get_ln_erk_timestepper() const
  {
    return timestepper_ln_erk;
  }

  // Getter for the simulationVariables object
  SimulationVariables* get_simulation_variables() const 
  { 
    return simVars;
  }

  // Getter for the number of levels
  int get_number_of_levels() const 
  {
    return levelSingletons->size();
  }
	      
  // Save the physical invariants
  void save_physical_invariants(
				int i_niter
				) 
  {
    time.push_back(simVars->timecontrol.current_timestep_size * i_niter);
    mass.push_back(simVars->diag.total_mass);
    energy.push_back(simVars->diag.total_energy);
    potentialEnstrophy.push_back(simVars->diag.total_potential_enstrophy);
  }

  // Getters for the time and invariants vectors
  const std::vector<double>& get_time()                const { return time; }
  const std::vector<double>& get_mass()                const { return mass; }
  const std::vector<double>& get_energy()              const { return energy; }
  const std::vector<double>& get_potential_enstrophy() const { return potentialEnstrophy; }

  // Getters for the residuals
  const std::vector<std::vector<double> >& get_residuals() const { return residuals; }
  std::vector<std::vector<double> >&       get_residuals()       { return residuals; }
    
protected:

  // Pointer to the SimulationVariables object
  SimulationVariables *simVars;

  // Pointer to the LevelSingletons vector
  std::vector<LevelSingleton> *levelSingletons;

  SWE_Sphere_TS_ln_erk*          timestepper_ln_erk;

  // Saved Residuals for each processor
  std::vector<std::vector<double> > residuals;

  // Diagnostics (mass, energy, enstrophy)
  SphereHelpers_Diagnostics* sphereDiagnostics;

  // Some contructors and operator= are disabled
  SphereDataCtx() {};
  SphereDataCtx(const SphereDataCtx&);
  SphereDataCtx& operator=(const SphereDataCtx&);
 
  // Vectors used for plotting
  std::vector<double> time;
  std::vector<double> mass;
  std::vector<double> energy;
  std::vector<double> potentialEnstrophy;

};

#endif
