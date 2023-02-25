/*
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_SHACKS_SHACK_PROG_ARG_DICT_HPP_
#define SRC_INCLUDE_SWEET_SHACKS_SHACK_PROG_ARG_DICT_HPP_


#include <list>
#include <memory>
#include <typeinfo>
#include <sweet/core/ErrorBase.hpp>
#include <sweet/core/ProgramArguments.hpp>
#include <sweet/core/shacks/ShackInterface.hpp>
#include <sweet/core/shacks/ShackDictionary.hpp>

namespace sweet
{

/*
 * A dictionary using class types as key.
 *
 * It also integrates parsing of program arguments.
 */
class ShackProgArgDictionary	:
		public ShackDictionary
{
private:
	ProgramArguments _programArguments;

	int _argc;
	char *const *_argv;

public:
	ShackProgArgDictionary(
			int i_argc,
			char *const *i_argv
	)	:
		ShackDictionary(),
		_argc(i_argc),
		_argv(i_argv)
	{
	}

	void setup()
	{
		_programArguments.setup(_argc, _argv);
		ERROR_FORWARD(_programArguments);
	}

public:
	void clear()
	{
		_programArguments.clear();
		ShackDictionary::clear();
	}

public:
	bool processProgramArguments()
	{
		return ShackDictionary::processProgramArguments(_programArguments);
	}


	bool processHelpArguments(bool i_with_error = true)
	{
		/*
		 * First, check for --help or -h
		 */
		if (_programArguments.argumentWithKeyExists("-h") || _programArguments.argumentWithKeyExists("--help"))
		{
			std::cout << "Printing help:" << std::endl;
			printProgramArguments();
			if (i_with_error)
				error.set("Help requested");

			return false;
		}

		return true;
	}

	bool checkAllArgumentsProcessed(bool i_create_error = true)
	{
		_programArguments.checkAllArgumentsProcessed(i_create_error);
		ERROR_FORWARD_WITH_RETURN_BOOLEAN(_programArguments);
	}

public:
	~ShackProgArgDictionary()
	{
	}
};

}

#endif
