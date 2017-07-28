#ifndef _PLANE_DATA_CTX_HPP_
#define _PLANE_DATA_CTX_HPP_

#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"

#include "SWE_Plane_TS_l_irk_n_erk.hpp"
#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_ln_erk.hpp"
#include "SWE_Plane_TS_l_direct.hpp"
#include "SWE_Plane_TS_l_rexi_n_erk.hpp"

// Class containing the context necessary to evaluate the right-hand sides
// Currently only contains a pointer to the level singletons and the SimulationVariables object

class PlaneDataCtx {

public:
  
  // Contructor
  PlaneDataCtx(
               SimulationVariables *i_simVars,
	       std::vector<LevelSingleton> *i_singletons
	       ) 
    : simVars(i_simVars),
      levelSingletons(i_singletons)
  {
    if (!simVars) 
      FatalError("PlaneDataCtx: simVars pointer is NULL!");

    if (!levelSingletons) 
      FatalError("PlaneDataCtx: levelSingletons pointer is NULL!");
    
    // initialize the time steppers from SWEET
    if (simVars->libpfasst.use_rexi)
      {
	timestepper_l_rexi_n_erk.resize(levelSingletons->size());	
	timestepper_l_irk_n_erk.resize(0);

      }
    else 
      {
	timestepper_l_irk_n_erk.resize(levelSingletons->size());
	timestepper_l_rexi_n_erk.resize(0);

      }
    timestepper_l_erk.resize(levelSingletons->size());
    ref_timestepper.resize(levelSingletons->size());

    for (int level = 0; level < timestepper_l_irk_n_erk.size(); ++level) 
      {
	// select first order integration in time for explicit
	// and first order integration for implicit (only order currently supported)
	simVars->disc.timestepping_order  = 1; 
	simVars->disc.timestepping_order2 = 1; 
	simVars->disc.use_staggering      = false; 
		
	// these timesteppers contain the functions called by LibPFASST 
	if (simVars->libpfasst.use_rexi)
	  {
	    timestepper_l_rexi_n_erk[level] = 
	      new SWE_Plane_TS_l_rexi_n_erk(
					   *simVars,
					   ((*levelSingletons)[level].op)
					   );
	    timestepper_l_rexi_n_erk[level]->get_implicit_timestepper().setup(simVars->rexi);
	  }
	else 
	  {
	    timestepper_l_irk_n_erk[level] = 
	      new SWE_Plane_TS_l_irk_n_erk(
					   *simVars,
					   ((*levelSingletons)[level].op)
					   );
	  }
     
	timestepper_l_erk[level] = 
	  new SWE_Plane_TS_l_erk(
				 *simVars,
				 ((*levelSingletons)[level].op)
				 );

	// use fourth order integration for the linear - nonlinear erk
	// this timestepper will be used to obtain a SWEET-generated reference solution
	simVars->disc.timestepping_order  = 4; 
	simVars->disc.timestepping_order2 = 4; 
	ref_timestepper[level] = 
	  new SWE_Plane_TS_ln_erk(
				  *simVars,
				  ((*levelSingletons)[level].op)
				  );
      }

  }

  // Destructor
  ~PlaneDataCtx() 
  {
    for (int level = 0; level < timestepper_l_irk_n_erk.size(); ++level) 
    {
      if (simVars->libpfasst.use_rexi) 
	delete timestepper_l_rexi_n_erk[level];
      else
	delete timestepper_l_irk_n_erk[level];
	
      delete timestepper_l_erk[level];
      delete ref_timestepper[level];
    }
  }

  // Getter for the plane data configuration at level i_level
  PlaneDataConfig* get_plane_data_config(
					 int i_level
					 ) const 
  {
    return &((*levelSingletons)[i_level].dataConfig);
  }

  // Getter for the plane data operators at level i_level
  PlaneOperators* get_plane_operators(
				      int i_level
				      ) const
  {
    return &((*levelSingletons)[i_level].op);
  }

  // Getter for the linear implicit nonlinear explicit SWEET time stepper at level i_level
  SWE_Plane_TS_l_irk_n_erk* get_l_irk_n_erk_timestepper(
							int i_level
							) const
  {
    if (simVars->libpfasst.use_rexi)
      return NULL;
    else 
      return timestepper_l_irk_n_erk[i_level];
  }

    // Getter for the linear implicit nonlinear explicit SWEET time stepper at level i_level
  SWE_Plane_TS_l_rexi_n_erk* get_l_rexi_n_erk_timestepper(
							int i_level
							) const
  {
    if (!simVars->libpfasst.use_rexi)
      return NULL;
    else 
      return timestepper_l_rexi_n_erk[i_level];
  }


  // Getter for the linear explicit SWEET time stepper at level i_level
  SWE_Plane_TS_l_erk* get_l_erk_timestepper(
						  int i_level
						  ) const
  {
    return timestepper_l_erk[i_level];
  }

  // Getter for the reference SWEET time stepper at level i_level
  SWE_Plane_TS_ln_erk* get_reference_timestepper(
						 int i_level
						 ) const
  {
    return ref_timestepper[i_level];
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

  // Pointer to the SWE_Plane time integrator (implicit linear part, explicit nonlinear part)
  std::vector<SWE_Plane_TS_l_irk_n_erk*>  timestepper_l_irk_n_erk;
  std::vector<SWE_Plane_TS_l_rexi_n_erk*> timestepper_l_rexi_n_erk;
  std::vector<SWE_Plane_TS_l_erk*>        timestepper_l_erk;

  // Pointer to the SWE_Plane time integrator used to compute a reference solution
  // this is currently ln_erk but might change later 
  std::vector<SWE_Plane_TS_ln_erk*> ref_timestepper;

  // Some contructors and operator= are disabled
  PlaneDataCtx() {};
  PlaneDataCtx(const PlaneDataCtx&);
  PlaneDataCtx& operator=(const PlaneDataCtx&);
 
};

#endif
