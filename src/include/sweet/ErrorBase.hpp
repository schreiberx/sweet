/*
 * ErrorBase.hpp
 *
 *  Created on: Feb 18, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_ERRORBASE_HPP_
#define SRC_INCLUDE_SWEET_ERRORBASE_HPP_

#include <string>
#include <ostream>
#include <sweet/SWEETError.hpp>


// For ptrace
#include <stdio.h>
#include <execinfo.h>

#include <sweet/Backtrace.hpp>

namespace sweet
{

class ErrorBase
{
	bool _hasError;
	std::string _errorMessage;


private:
	/*
	 * Check for environment variable SWEET_ERROR_WITH_STACKTRACE and if it exists,
	 * add stack trace to error information
	 */
	bool _withStacktrace()
	{
		static bool envLoaded = false;
		static bool envExists = false;

		if (!envLoaded)
		{
			char *val = std::getenv("SWEET_ERROR_WITH_STACKTRACE");
			envExists = (val != nullptr);
			envLoaded = true;
		}
		return envExists;
	}


public:
	ErrorBase()	:
		_hasError(false)
	{
	}

	/*
	 * Set an error
	 *
	 * \return always *false* to be able to use a one-liner error message
	 */
	bool set(const std::string &i_errorMessage)
	{
		if (_hasError)
		{
			_errorMessage += " | " + i_errorMessage;
		}
		else
		{
			_hasError = true;
			_errorMessage = i_errorMessage;
		}


		if (_withStacktrace())
		{
			std::string gdb = Backtrace::getGDBBacktrace();

			if (gdb != "")
			{
				_errorMessage += "\n";
				_errorMessage += "Stacktrace (from GDB):\n";
				_errorMessage += Backtrace::getGDBBacktrace();
			}
		}
		return false;
	}

	bool forwardFrom(ErrorBase &i_error)
	{
		if (!i_error._hasError)
			return false;

		_hasError = i_error._hasError;
		_errorMessage = i_error._errorMessage;

		i_error.reset();

		return true;
	}

	/*
	 * Forward errors
	 *
	 * \return **true** if there's no error
	 */
	bool forwardFromWithPositiveReturn(ErrorBase &i_error)
	{
		if (!i_error._hasError)
			return true;

		_hasError = i_error._hasError;
		_errorMessage = i_error._errorMessage;

		i_error.reset();

		return false;
	}

	bool exists() const
	{
		return _hasError;
	}

	void assertNoError() const
	{
		if (_hasError)
		{
			std::cerr << "Internal problem detected" << std::endl;
			SWEETError("Use backtrace to debug this problem!");
			exit(1);
		}
	}

	void reset()
	{
		_hasError = false;
		_errorMessage = "";
	}

	std::string get()
	{
		if (!_hasError)
		{
			std::cerr << "Error message requested, but has no error!" << std::endl;
			SWEETError("Use backtrace to debug this problem!");
			exit(1);
		}

		std::string tmp = _errorMessage;
		reset();
		return tmp;
	}

public:
	void print(std::ostream &io_os = std::cerr)
	{
		io_os << "ERROR: " << get() << std::endl;
	}

	~ErrorBase()
	{
		if (_hasError)
		{
			std::cerr << "ERROR was not processed! Error message: '" << _errorMessage << "'" << std::endl;
			exit(1);
		}
	}
};

}

#endif
