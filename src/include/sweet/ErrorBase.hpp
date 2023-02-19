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

	bool _errorProcessed;

public:
	ErrorBase()	:
		_hasError(false),
		_errorProcessed(true)
	{
	}

	void errorSet(const std::string &i_errorMessage)
	{
		_hasError = true;
		_errorMessage = i_errorMessage;
	}

	bool errorForward(ErrorBase &i_errorBase)
	{
		if (!i_errorBase._hasError)
			return false;

		_hasError = i_errorBase._hasError;
		_errorMessage = i_errorBase._errorMessage;
		_errorProcessed = i_errorBase._errorProcessed;

		i_errorBase.errorReset();

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
		_errorProcessed = true;
	}

	const std::string& errorGet()
	{
		if (!_hasError)
		{
			std::cerr << "Error message requested, but has no error!" << std::endl;
			exit(1);
		}
		_errorProcessed = true;
		return _errorMessage;
	}

	~ErrorBase()
	{
		if (_hasError && !_errorProcessed)
		{
			std::cerr << "ERROR: Lack of processing error '" << _errorMessage << "'" << std::endl;
			exit(1);
		}
	}
};

}

#endif
