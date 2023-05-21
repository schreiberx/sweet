/*
 * ClassDictionaryInterface.hpp
 *
 *  Created on: Feb 19, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKINTERFACE_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKINTERFACE_HPP_


#include <typeinfo>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/ErrorBase.hpp>
#include <string>

namespace sweet
{

class ShackInterface
{
public:
	ErrorBase error;

	bool argumentsProcessed;

	ShackInterface()	:
		argumentsProcessed(false)
	{
	}

	virtual void printProgramArguments(const std::string &i_prefix = "")
	= 0;
	//{}

	virtual bool processProgramArguments(ProgramArguments &i_pa)
	= 0;
	//{	return true;}

	virtual void printShack(const std::string &i_prefix = "")
	= 0;
	//{	std::cout << "[TODO]" << std::endl; }


	virtual ~ShackInterface()
	{}
};

}

#endif
