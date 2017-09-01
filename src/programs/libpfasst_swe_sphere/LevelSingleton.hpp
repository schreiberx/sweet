#ifndef _LEVEL_SINGLETON_CTX_HPP_
#define _LEVEL_SINGLETON_CTX_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereOperators.hpp>

// Class to store the configurations and operators at each level

class LevelSingleton
{

public:
  
  int level;
  SphereDataConfig dataConfig;
  SphereOperators op;
};

#endif
