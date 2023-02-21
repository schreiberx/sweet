/*
 * ShackSDC.hpp
 *
 *  Created on: Feb 21, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKSDC_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKSDC_HPP_

#include <string>
#include <iostream>
#include <sweet/ProgramArguments.hpp>
#include <sweet/shacks/ShackInterface.hpp>


/**
 * SDC parameters
 */
class SDC	:
		public sweet::ClassDictionaryInterface
{
public:
	std::string fileName;

	// Nodes values (between 0 and 1)
	sweet::DictArrayND<1, double> nodes;
	size_t nNodes=0;

	// Quadrature weights
	sweet::DictArrayND<1, double> weights;

	// Collocation matrix
	sweet::DictArrayND<2, double> qMatrix;

	// QDelta matrix for implicit sweep
	sweet::DictArrayND<2, double> qDeltaI;

	// QDelta matrix for explicit sweep
	sweet::DictArrayND<2, double> qDeltaE;

	// QDelta matrix for initial (implicit) sweep
	sweet::DictArrayND<2, double> qDelta0;

	// Number of iterations (sweeps)
	sweet::Dict::int64 nIter=0;

	// Type of initial sweep to use
	std::string initSweepType="COPY";

	// Wether or not use the diagonal implementation
	sweet::Dict::int64 diagonal=0;

	// Wether or not use collocation update for end point
	sweet::Dict::int64 useEndUpdate=0;


	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << "SDC option:" << std::endl;
		std::cout << "	--sdc-file [path]   SDC parameters in sweet::Dict format" << std::endl;
		std::cout << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		if (i_pa.getArgumentValueByKey("--sdc-file", fileName))
		{
			sweet::Dict params(fileName);
			params.get("nodes", nodes);
			nNodes = nodes.size();
			params.get("weights", weights);
			params.get("qMatrix", qMatrix);
			params.get("qDeltaI", qDeltaI);
			params.get("qDeltaE", qDeltaE);
			params.get("qDelta0", qDelta0);
			params.get("nIter", nIter);
			params.get("diagonal", diagonal);
			params.get("initSweepType", initSweepType);
			params.get("useEndUpdate", useEndUpdate);
		}

		return error.forwardFromWithPositiveReturn(i_pa.error);
	}

	virtual void printClass(
		const std::string& i_prefix = ""
	)
	{
		std::cout << std::endl;
		std::cout << "SDC:" << std::endl;
		std::cout << " + M (number of nodes): " << nodes.size() << std::endl;
		std::cout << " + nodes: " << nodes << std::endl;
		std::cout << " + weights: " << weights << std::endl;
		std::cout << " + qMatrix: " << qMatrix << std::endl;
		std::cout << " + qDeltaI: " << qDeltaI << std::endl;
		std::cout << " + qDeltaE: " << qDeltaE << std::endl;
		std::cout << " + qDelta0: " << qDelta0 << std::endl;
		std::cout << " + nIter: " << nIter << std::endl;
		std::cout << " + diagonal: " << diagonal << std::endl;
		std::cout << " + qDeltaInit: " << diagonal << std::endl;
		std::cout << " + useEndUpdate: " << diagonal << std::endl;
		std::cout << std::endl;
	}
};





#endif /* SRC_INCLUDE_SWEET_SHACKS_SHACKSDC_HPP_ */
