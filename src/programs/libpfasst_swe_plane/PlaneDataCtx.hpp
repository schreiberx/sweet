#ifndef _PLANE_DATA_CTX_HPP_
#define _PLANE_DATA_CTX_HPP_

#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"
#include "SWE_Plane_TS_l_erk.hpp"
#include "SWE_Plane_TS_l_direct.hpp"

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

    // select fourth order integration in time
    simVars->disc.timestepping_order = 4; 

    // initialize the time stepper from SWEET
    timestepping_interfaces.resize(levelSingletons->size());
    ref_timestepping_interfaces.resize(levelSingletons->size());

    for (int level = 0; level < timestepping_interfaces.size(); ++level) 
    {
      timestepping_interfaces[level] = 
	new SWE_Plane_TS_l_erk(
	                       *simVars,
			       ((*levelSingletons)[level].op)
		              );
      ref_timestepping_interfaces[level] = 
	new SWE_Plane_TS_l_direct(
				  *simVars,
				  ((*levelSingletons)[level].op)
				  );

    }
  }

  // Destructor
  ~PlaneDataCtx() 
  {
    for (int level = 0; level < timestepping_interfaces.size(); ++level) 
    {
      delete timestepping_interfaces[level];
      delete ref_timestepping_interfaces[level];
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

  // Getter for the SWEET time stepper at level i_level
  SWE_Plane_TS_interface* get_timestepper(
                                          int i_level
					  ) const
  {
    return timestepping_interfaces[i_level];
  }

  // Getter for the reference SWEET time stepper at level i_level
  SWE_Plane_TS_interface* get_reference_timestepper(
                                                    int i_level
					           ) const
  {
    return ref_timestepping_interfaces[i_level];
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

  // Pointer to the SWE_Plane time integrator
  std::vector<SWE_Plane_TS_interface*> timestepping_interfaces;

  // Pointer to the SWE_Plane time integrator used to compute a reference solution
  std::vector<SWE_Plane_TS_interface*> ref_timestepping_interfaces;

  // Some contructors and operator= are disabled
  PlaneDataCtx() {};
  PlaneDataCtx(const PlaneDataCtx&);
  PlaneDataCtx& operator=(const PlaneDataCtx&);
 
};

#endif
