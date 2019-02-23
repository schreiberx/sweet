#ifndef _LEVEL_SINGLETON_CTX_HPP_
#define _LEVEL_SINGLETON_CTX_HPP_

#include <sweet/sphere/SphereDataSpectral.hpp>
#include <sweet/sphere/SphereOperators.hpp>
#include <benchmarks_sphere/SWESphereBenchmarksCombined.hpp>

// Class to store the configurations and operators at each level

class LevelSingleton
{

public:
  
  int              level;
  SphereDataConfig dataConfig;
  SphereDataConfig dataConfigNoDealiasing;
  SphereOperators  op;
  SphereOperators  opNoDealiasing;

  SWESphereBenchmarksCombined benchmarks;
};

#endif
