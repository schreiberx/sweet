#ifndef _SPHERE_DATA_CTX_HPP_
#define _SPHERE_DATA_CTX_HPP_

#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"

#include "SWE_Sphere_TS_l_erk.hpp"
#include "SWE_Sphere_TS_l_erk_n_erk.hpp"
#include "SWE_Sphere_TS_l_irk.hpp"
#include "SWE_Sphere_TS_l_rexi.hpp"

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
    if (!simVars) 
      FatalError("SphereDataCtx: simVars pointer is NULL!");

    if (!levelSingletons) 
      FatalError("SphereDataCtx: levelSingletons pointer is NULL!");
    
    // initialize the time steppers from SWEET
    timestepper_l_erk_n_erk.resize(levelSingletons->size());
    if (simVars->libpfasst.use_rexi) 
      {
	timestepper_l_irk.resize(levelSingletons->size());
      }
    else
      {
	timestepper_l_rexi.resize(levelSingletons->size());
      }

    for (int level = 0; level < levelSingletons->size(); ++level) 
      {
	// select first order integration in time for explicit
	// and first order integration for implicit (only order currently supported)
	simVars->disc.timestepping_order  = 1; 
	simVars->disc.timestepping_order2 = 1; 
		
	// these timesteppers contain the functions called by LibPFASST 
	timestepper_l_erk_n_erk[level] = 
	  new SWE_Sphere_TS_l_erk_n_erk(
				   *simVars,
				   ((*levelSingletons)[level].op)
				   );
        timestepper_l_erk_n_erk[level]->setup(simVars->disc.timestepping_order);

        if (simVars->libpfasst.use_rexi) 
	  {
            timestepper_l_rexi[level] = 
	      new SWE_Sphere_TS_l_rexi(
	                               *simVars,
				       ((*levelSingletons)[level].op)
				       );

            timestepper_l_rexi[level]->setup(simVars->rexi,
					     "phi0",
					     simVars->timecontrol.current_timestep_size,
					     simVars->sim.f_sphere,
					     (simVars->sim.coriolis_omega == 0)
					     );
	  }
	else
	  {
	    timestepper_l_irk[level] = 
	      new SWE_Sphere_TS_l_irk(
				      *simVars,
				      ((*levelSingletons)[level].op)
				      );
	
	    timestepper_l_irk[level]->setup(simVars->disc.timestepping_order,
					    simVars->timecontrol.current_timestep_size,
					    simVars->rexi.use_sphere_extended_modes);
	  }
      }
  }

  // Destructor
  ~SphereDataCtx() 
  {
    for (int level = 0; level < timestepper_l_erk_n_erk.size(); ++level) 
    {
      if (simVars->libpfasst.use_rexi)
	{
	  delete timestepper_l_rexi[level];
	}
      else 
	{
	  delete timestepper_l_irk[level];
	}
      delete timestepper_l_erk_n_erk[level];
    }
  }

  // Getter for the sphere data configuration at level i_level
  SphereDataConfig* get_sphere_data_config(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].dataConfig);
  }

  // Getter for the sphere data operators at level i_level
  SphereOperators* get_sphere_operators(
				      int i_level
				      ) const
  {
    return &((*levelSingletons)[i_level].op);
  }

  // Getter for the linear implicit nonlinear explicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_erk_n_erk* get_l_erk_n_erk_timestepper(
					       int i_level
					       ) const
  {
    return timestepper_l_erk_n_erk[i_level];
  }

  // Getter for the linear explicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_irk* get_l_irk_timestepper(
					     int i_level
					     ) const
  {
    if (simVars->libpfasst.use_rexi)
      return NULL;
    else
      return timestepper_l_irk[i_level];
  }

  // Getter for the REXI implicit SWEET time stepper at level i_level
  SWE_Sphere_TS_l_rexi* get_l_rexi_timestepper( 
					       int i_level
					       ) const
  {
    if (!simVars->libpfasst.use_rexi)
      return NULL;
    else
      return timestepper_l_rexi[i_level];
  }
 
  // Getter for the simulationVariables object
  SimulationVariables* get_simulation_variables() const 
  { 
    return simVars;
  }

  // Setter for the timestepping params
  void setup_time_steps(
			double i_t, 
			double i_dt
			)
  {
    simVars->timecontrol.current_timestep_size   = i_dt;
    simVars->timecontrol.current_simulation_time = i_t;
    simVars->timecontrol.max_simulation_time     = i_dt;
  }
	      
    
protected:

  // Pointer to the SimulationVariables object
  SimulationVariables *simVars;

  // Pointer to the LevelSingletons vector
  std::vector<LevelSingleton> *levelSingletons;

  // Pointer to the SWE_Sphere time integrator (implicit linear part, explicit nonlinear part)
  std::vector<SWE_Sphere_TS_l_erk_n_erk*> timestepper_l_erk_n_erk;
  std::vector<SWE_Sphere_TS_l_irk*> timestepper_l_irk;
  std::vector<SWE_Sphere_TS_l_rexi*> timestepper_l_rexi; 

  // Some contructors and operator= are disabled
  SphereDataCtx() {};
  SphereDataCtx(const SphereDataCtx&);
  SphereDataCtx& operator=(const SphereDataCtx&);
 
};

#endif
