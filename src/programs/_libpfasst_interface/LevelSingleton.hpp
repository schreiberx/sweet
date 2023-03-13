#ifndef _LEVEL_SINGLETON_CTX_HPP_
#define _LEVEL_SINGLETON_CTX_HPP_

#include "../pde_sweSphere/PDESWESphere_BenchmarksCombined.hpp"
#include <sweet/core/sphere/SphereData_Spectral.hpp>
#include <sweet/core/sphere/SphereOperators.hpp>

// Class to store the configurations and operators at each level

class LevelSingleton
{

public:
  
  int              level;
  sweet::SphereData_Config dataConfig;
  sweet::SphereOperators  op;

  PDESWESphere_BenchmarksCombined benchmarks;
};

#endif
