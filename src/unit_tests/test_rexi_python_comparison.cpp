/*
 * test_rexi.cpp
 *
 *  Created on: 2 Aug 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#include <iostream>
#include <rexi/REXI.hpp>
#include <sweet/SimulationVariables.hpp>


#if 0
	typedef __float128 TGeneration;
#else
	typedef double TGeneration;
#endif
typedef std::complex<TGeneration> TComplexGeneration;



typedef double T;
typedef std::complex<T> TComplex;



int main(
		int i_argc,
		char *const i_argv[]
)
{
	SimulationVariables simVars;
	if (!simVars.setupFromMainParameters(i_argc, i_argv, nullptr, false))
	{
		return -1;
	}

	REXI_Terry<TGeneration, T> rexi("phi0", 0.2, 2, 0, false, false);
	for (int i = 0; i < rexi.beta_re.size(); i++)
	{
		std::cout << rexi.beta_re[i] << std::endl;
	}

	return 0;
}
