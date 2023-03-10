#ifndef _LEVEL_SINGLETON_CTX_HPP_
#define _LEVEL_SINGLETON_CTX_HPP_

#include "../swe_sphere_benchmarks/BenchmarksSphereSWE.hpp"
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>

// Class to store the configurations and operators at each level

class LevelSingleton
{

public:
  
  int              level;
  sweet::SphereData_Config dataConfig;
  SphereOperators  op;

  BenchmarksSphereSWE benchmarks;
};

#endif
