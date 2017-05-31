#ifndef _PLANE_DATA_CTX_HPP_
#define _PLANE_DATA_CTX_HPP_

#include <vector>
#include <sweet/SimulationVariables.hpp>
#include "LevelSingleton.hpp"

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

  // Getter for the simulationVariables object
  SimulationVariables* get_simulation_variables() const 
  { 
    return simVars;
  }
    
protected:

  // Pointer to the SimulationVariables object
  SimulationVariables *simVars;

  // Pointer to the LevelSingletons vector
  std::vector<LevelSingleton> *levelSingletons;

  // Some contructors and operator= are disabled
  PlaneDataCtx() {};
  PlaneDataCtx(const PlaneDataCtx&);
  PlaneDataCtx& operator=(const PlaneDataCtx&);
 
};

#endif
