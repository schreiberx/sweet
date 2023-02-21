/*
 * ClassDictionaryInterface.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com> Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACKINTERFACE_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACKINTERFACE_HPP_


#include <typeinfo>
#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include <string>

namespace sweet
{

class ClassDictionaryInterface
{
public:
	sweet::ErrorBase error;

	virtual void printProgramArguments(const std::string &i_prefix = "")
	= 0;
	//{}

	virtual bool processProgramArguments(sweet::ProgramArguments &i_pa)
	= 0;
	//{	return true;}

	virtual void printClass(const std::string &i_prefix = "")
	= 0;
	//{	std::cout << "[TODO]" << std::endl; }


	virtual ~ClassDictionaryInterface()
	{}
};

}

#endif
