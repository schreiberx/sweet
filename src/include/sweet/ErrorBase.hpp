/*
 * ErrorBase.hpp
 *
 *  Created on: Feb 18, 2023
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_ERRORBASE_HPP_
#define SRC_INCLUDE_SWEET_ERRORBASE_HPP_

#include <string>
#include <iostream>

namespace sweet
{

class ErrorBase
{
	bool _hasError;
	std::string _errorMessage;

public:
	ErrorBase()	:
		_hasError(false)
	{
	}

	void errorSet(const std::string &i_errorMessage)
	{
		_hasError = true;
		_errorMessage = i_errorMessage;
	}

	bool errorForwardFrom(ErrorBase &i_error)
	{
		if (!i_error._hasError)
			return false;

		_hasError = i_error._hasError;
		_errorMessage = i_error._errorMessage;

		i_error.errorReset();

		return true;
	}

	bool errorExists() const
	{
		return _hasError;
	}

	void errorAssertNoError() const
	{
		if (_hasError)
		{
			std::cout << "Internal problem detected" << std::endl;
			exit(1);
		}
	}

	void errorReset()
	{
		_hasError = false;
		_errorMessage = "";
	}

	std::string errorGet()
	{
		if (!_hasError)
		{
			std::cerr << "Error message requested, but has no error!" << std::endl;
			exit(1);
		}

		std::string tmp = _errorMessage;
		errorReset();
		return tmp;
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
