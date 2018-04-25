#ifndef _LEVEL_SINGLETON_CTX_HPP_
#define _LEVEL_SINGLETON_CTX_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneOperators.hpp>

// Class to store the configurations and operators at each level

class LevelSingleton
{

public:
  
  int             level;
  PlaneDataConfig dataConfig;
  PlaneDataConfig dataConfigNoDealiasing;
  PlaneOperators  op;
  PlaneOperators  opNoDealiasing;
};

#endif
