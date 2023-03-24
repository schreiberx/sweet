/*
 * ShackSDC.hpp
 *
 *  Created on: Feb 21, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKSDC_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKSDC_HPP_

#include <string>
#include <iostream>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/dict/Dict.hpp>

namespace sweet {

class ShackSDC	:
		public sweet::ShackInterface
{
public:
	std::string fileName;

	// Nodes values (between 0 and 1)
	sweet::DictArrayND<1, double> nodes;

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
	sweet::Dict::int64 nIter;

	// Type of initial sweep to use
	std::string initialSweepType;

	// Wether or not use the diagonal implementation
	sweet::Dict::int64 diagonal;

	// Wether or not activate time parallelization
	int parallel;

	// Wether or not use collocation update for end point
	sweet::Dict::int64 useEndUpdate;

	// Unique string ID
	std::string idString;

	
	// Default parameters for SDC shack
	ShackSDC(){
		int nNodes = 3;
		nIter = 3;
		diagonal = false;
		parallel = false;
		initialSweepType = "COPY";
		useEndUpdate = false;
		idString = "M3_RADAU-RIGHT_K3_BE_FE_COPY";

		// RADAU-RIGHT nodes, weights quadrature matrix
		nodes.setup(nNodes);
		const double _nodes[] = {
			0.15505102572168, 0.64494897427832, 1.
		};
		nodes = _nodes;

		weights.setup(nNodes);
		const double _weights[] = {
			0.3764030627004656, 0.51248582618842650, 0.1111111111111111
		};
		weights = _weights;

		qMatrix.setup(nNodes, nNodes);
		const double _qMatrix[] = {
			0.1968154772236567, -0.06553542585019642,  0.02377097434821968,
			0.394424314739085,   0.2920734116652353,  -0.04154875212600038,
			0.3764030627004656,  0.5124858261884265,   0.1111111111111079
		};
		qMatrix = _qMatrix;

		// BE for implicit sweep
		qDeltaI.setup(nNodes, nNodes);
		const double _qDeltaI[] = {
			0.15505102572168, 0.,         0.,
 			0.15505102572168, 0.48989794855664, 0.,
			0.15505102572168, 0.48989794855664, 0.35505102572168        
		};
		qDeltaI = _qDeltaI;

		// FE for explicit sweep
		qDeltaE.setup(nNodes, nNodes);
		const double _qDeltaE[] = {
			0.,         0.,         0.,
 			0.48989794855664, 0.,         0.,
			0.48989794855664, 0.35505102572168, 0.        
		};
		qDeltaE = _qDeltaE;

		// BEpar for initial sweep
		qDelta0.setup(nNodes, nNodes);
		const double _qDelta0[] = {
			0.15505102572168, 0.		      , 0.,
 			0.		        , 0.64494897427832, 0.,
			0.		        , 0.		      , 1.        
		};
		qDelta0 = _qDelta0;

	}

	void printProgramArguments(const std::string& i_prefix = "")
	{
		std::cout << std::endl;
		std::cout << "SDC option:" << std::endl;
		std::cout << "	--sdc-file 	   [path]  SDC parameters in Dict format" << std::endl;
		std::cout << "	--sdc-parallel [bool]  wether or not activate time parallelization" << std::endl;
		std::cout << std::endl;
	}

	bool processProgramArguments(sweet::ProgramArguments &i_pa)
	{
		i_pa.getArgumentValueByKey("--sdc-parallel", parallel);

		if (i_pa.getArgumentValueByKey("--sdc-file", fileName))
		{
			sweet::Dict params(fileName);
			params.get("nodes", nodes);
			params.get("weights", weights);
			params.get("qMatrix", qMatrix);
			params.get("qDeltaI", qDeltaI);
			params.get("qDeltaE", qDeltaE);
			params.get("qDelta0", qDelta0);
			params.get("nIter", nIter);
			params.get("diagonal", diagonal);
			params.get("initialSweepType", initialSweepType);
			params.get("useEndUpdate", useEndUpdate);
			params.get("idString", idString);
		}

		return error.forwardWithPositiveReturn(i_pa.error);
	}

	virtual void printShack(
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
		std::cout << " + parallel: " << parallel << std::endl;
		std::cout << " + initialSweepType: " << initialSweepType << std::endl;
		std::cout << " + useEndUpdate: " << useEndUpdate << std::endl;
		std::cout << std::endl;
	}
};


}


#endif
