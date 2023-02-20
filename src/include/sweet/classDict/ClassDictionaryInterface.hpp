/*
 * ClassDictionaryInterface.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_CLASSDICT_CLASSDICTIONARYINTERFACE_HPP_
#define SRC_INCLUDE_SWEET_CLASSDICT_CLASSDICTIONARYINTERFACE_HPP_


#include <typeinfo>
#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include <string>

class ClassDictionaryInterface
{
public:
	virtual sweet::ErrorBase& getError() = 0;

	virtual void outputProgramArguments(std::string i_prefix = "") = 0;

	virtual bool processProgramArguments(sweet::ProgramArguments &i_pa) = 0;

	virtual void outputVariables(std::string i_prefix = "") = 0;


	virtual ~ClassDictionaryInterface()
	{

	}
};


#endif
