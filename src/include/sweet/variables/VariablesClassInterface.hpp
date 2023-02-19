/*
 * VariablesClassInterface.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_VARIABLES_VARIABLESCLASSINTERFACE_HPP_
#define SRC_INCLUDE_SWEET_VARIABLES_VARIABLESCLASSINTERFACE_HPP_


#include <typeinfo>
#include <sweet/ProgramArguments.hpp>
#include <sweet/ErrorBase.hpp>
#include <string>

class VariablesClassDictionaryInterface
{
public:
	virtual sweet::ErrorBase& getError() = 0;

	virtual void outputProgramArguments(std::string i_prefix = "") = 0;

	virtual bool processProgramArguments(sweet::ProgramArguments &i_pa) = 0;

	virtual void outputVariables(std::string i_prefix = "") = 0;


	virtual ~VariablesClassDictionaryInterface(){}
};


#endif
