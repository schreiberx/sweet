/*
 * ProgramArguments.hpp
 *
 *  Created on: Feb 18, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_PROGRAMARGUMENTS_HPP_
#define SRC_INCLUDE_SWEET_PROGRAMARGUMENTS_HPP_

#include <vector>
#include <list>
#include <string>
#include <ostream>
#include <iostream>
#include <sstream>

#include <sweet/ErrorBase.hpp>

namespace sweet
{

/*
 * Program argument parser
 *
 * Get arguments of the form
 *
 * -T 123
 * --foo-bar=123
 *
 * Supported, but not recommended:
 * --foo-foo message
 *
 * The **parameters** are the values
 */
class ProgramArguments
{
public:
	ErrorBase error;

private:
	class _ProgramArgument
	{
	public:
		std::string key;
		std::string value;
		bool argumentParsedAndAccessed;

		_ProgramArgument(
				const std::string &i_key,
				const std::string &i_value
		)	:
			key(i_key),
			value(i_value),
			argumentParsedAndAccessed(false)
		{
		}


		_ProgramArgument()	:
			argumentParsedAndAccessed(false)
		{
		}

		void reset()
		{
			key = "";
			value = "";
			argumentParsedAndAccessed = false;
		}
	};

	std::string _arg0;
	std::vector<_ProgramArgument> _arguments;

	bool _errorForDoubleParsing;

	bool _errorForDuplicateKeysInParsing;

	std::list<std::string> _keysInParsing;

public:

	ProgramArguments(
			bool i_errorForDoubleParsing = true,	///< Trigger an error if an argument is parsed twice
			bool i_errorForDuplicateKeysInParsing = true		///< Trigger an error if there are duplicate keys
	)	:
		_errorForDoubleParsing(i_errorForDoubleParsing),
		_errorForDuplicateKeysInParsing(i_errorForDuplicateKeysInParsing)
	{
	}


private:
	void _addArguments(
			const std::string i_key,
			const std::string i_value
	)
	{
		if (argumentWithKeyExists(i_key))
		{
			error.errorSet("Argument with key '"+i_key+"' already exists - did you specify this twice?");
			return;
		}
		_arguments.push_back(
			_ProgramArgument(
					i_key,
					i_value
				)
			);
	}

public:
	bool setup(int argc, const char *const argv[])
	{
		_arg0 = argv[0];

		/*
		 * 0: no previous processed data
		 * 1: key processed
		 */
		int state = 0;

		std::string tmp_key;

		for (int i = 1; i < argc; i++)
		{
			std::string arg = argv[i];

			if (state == 0)
			{
				/*
				 * Check for single or double dash
				 */

				if (arg.size() <= 1)
				{
					std::stringstream ss;
					ss << "Error parsing argument " << i << ": '" << arg << "' (too short)" << std::endl;
					error.errorSet(ss.str());
					return false;
				}

				int num_dashes = 0;
				if (arg[0] == '-')
				{
					if (arg[1] == '-')
					{
						num_dashes = 2;
					}
					else
					{
						num_dashes = 1;
					}
				}

				if (num_dashes == 0)
				{
					std::stringstream ss;
					ss << "Error parsing argument " << i << ": '" << arg << "' (missing dashes)" << std::endl;
					error.errorSet(ss.str());
					return false;
				}

				/*
				 * Search for "=" separator
				 */
				std::size_t pos = arg.find('=');

				if (pos == std::string::npos)
				{
					/*
					 * Not found => continue
					 */
					state = 1;
					tmp_key = arg.substr(num_dashes, std::string::npos);
					continue;
				}

				_addArguments(
						arg.substr(num_dashes, pos-num_dashes),
						arg.substr(pos+1, std::string::npos)
					);
				continue;
			}

			if (state == 1)
			{
				_addArguments(tmp_key, arg);

				state = 0;
				continue;
			}
		}

		if (state == 1)
		{
			error.errorSet("Invalid format of program arguments (last one could not be parsed)");
			return false;
		}

		return true;
	}

public:
	bool argumentWithKeyExists(const std::string& i_key)
	{
		for (int i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (a.key == i_key)
			{
				return true;
			}
		}

		return false;
	}


public:
	bool checkAllArgumentsProcessed(bool i_create_error = true)
	{
		for (int i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (!a.argumentParsedAndAccessed)
			{
				if (i_create_error)
					error.errorSet("Argument '"+a.key+"' not processed!");

				return false;
			}
		}

		return true;
	}

	bool _getFullArgumentByKey(
			const std::string& i_key,
			_ProgramArgument** o_pa,
			bool i_no_error_if_key_is_missing = true
	)
	{
		for (int i = 0; i < _arguments.size(); i++)
		{
			_ProgramArgument &a = _arguments[i];
			if (a.key == i_key)
			{
				*o_pa = &a;
				return true;
			}
		}

		if (!i_no_error_if_key_is_missing)
			error.errorSet(std::string("")+"Key '"+i_key+"' not found");

		return false;
	}


	bool checkAndAddDuplicateKeys(const std::string& i_key)
	{
		for (std::list<std::string>::iterator i = _keysInParsing.begin(); i != _keysInParsing.end(); i++)
		//for (int i = 0; i < _keysInParsing.size(); i++)
		{
			if (i_key == *i)
			{
				error.errorSet("Key '"+i_key+"' parsed twice");
				return false;
			}
		}

		_keysInParsing.push_back(i_key);

		return true;
	}


	/**
	 * Get double precision value
	 */
public:
	bool getArgumentValueByKey(
			const std::string& i_key,
			double &o_value,
			bool i_no_error_if_key_is_missing = true
	)
	{
		if (_errorForDuplicateKeysInParsing)
			if (!checkAndAddDuplicateKeys(i_key))
				return false;

		_ProgramArgument *pa;
		if (!_getFullArgumentByKey(i_key, &pa))
		{
			error.errorReset();
			if (!_getFullArgumentByKey(i_key, &pa))
			{
				return false;
			}
		}

		if (_errorForDoubleParsing && pa->argumentParsedAndAccessed)
		{
			error.errorSet("Argument with key '"+i_key+"' already parsed");
			return false;
		}

		try
		{
			o_value = std::stod(pa->value);
		}
		catch (const std::exception &e)
		{
			error.errorSet("Exception caught during conversion of value '"+pa->value+"' to double: "+e.what());
			return false;
		}

		pa->argumentParsedAndAccessed = true;

		return true;
	}

	/**
	 * Get integer value
	 */
	bool getArgumentValueByKey(
			const std::string& i_key,
			int &o_value,
			bool i_no_error_if_key_is_missing = true
	)
	{
		_ProgramArgument *pa;
		if (!_getFullArgumentByKey(i_key, &pa))
		{
			error.errorReset();
			if (!_getFullArgumentByKey(i_key, &pa))
			{
				return false;
			}
		}

		if (_errorForDoubleParsing && pa->argumentParsedAndAccessed)
		{
			error.errorSet("Argument with key '"+i_key+"' already parsed");
			return false;
		}

		try
		{
			o_value = std::stoi(pa->value);
		}
		catch (const std::exception &e)
		{
			error.errorSet("Exception caught during conversion of value '"+pa->value+"' to double: "+e.what());
			return false;
		}

		pa->argumentParsedAndAccessed = true;

		return true;
	}

	/*
	 * Get boolean value
	 */
	bool getArgumentValueByKey(
			const std::string& i_key,
			bool &o_value,
			bool i_no_error_if_key_is_missing = true
	)
	{
		_ProgramArgument* pa;
		if (!_getFullArgumentByKey(i_key, &pa))
		{
			error.errorReset();
			if (!_getFullArgumentByKey(i_key, &pa))
			{
				return false;
			}
		}

		if (_errorForDoubleParsing && pa->argumentParsedAndAccessed)
		{
			error.errorSet("Argument with key '"+i_key+"' already parsed");
			return false;
		}

		std::string val = pa->value;

		for (int i = 0; i < val.length(); i++)
		{
			if (val[i] >= 'A' && val[i] <= 'Z')
				val[i] -= 'A'-'a';
		}

		if (val == "1" || val == "true")
		{
			o_value = true;
			pa->argumentParsedAndAccessed = true;
			return true;
		}

		if (val == "0" || val == "false")
		{
			o_value = false;
			pa->argumentParsedAndAccessed = true;
			return true;
		}

		error.errorSet("Cannot parse value '" + val +"' as boolean type");

		return false;
	}

	template <typename T>
	bool getArgumentValueBy2Keys(
			const std::string& i_key1,
			const std::string& i_key2,
			T &o_value,
			bool i_no_error_if_key_is_missing = true
	)
	{
		_ProgramArgument pa;

		error.errorAssertNoError();
		if (getArgumentValueByKey(i_key1, o_value, i_no_error_if_key_is_missing))
			return true;

		if (error.errorExists())
			return false;

		if (getArgumentValueByKey(i_key2, o_value, i_no_error_if_key_is_missing))
			return true;

		return false;
	}

	template <typename T>
	bool getArgumentValueBy3Keys(
			const std::string& i_key1,
			const std::string& i_key2,
			const std::string& i_key3,
			T &o_value,
			bool i_no_error_if_key_is_missing = true
	)
	{
		_ProgramArgument pa;

		error.errorAssertNoError();
		if (getArgumentValueByKey(i_key1, o_value, i_no_error_if_key_is_missing))
			return true;

		if (error.errorExists())
			return false;

		if (getArgumentValueByKey(i_key2, o_value, i_no_error_if_key_is_missing))
			return true;

		if (error.errorExists())
			return false;

		if (getArgumentValueByKey(i_key3, o_value, i_no_error_if_key_is_missing))
			return true;

		return false;
	}

	friend
	std::ostream&
	operator<<(std::ostream &io_os, const ProgramArguments i_pa)
	{
		for (int i = 0; i < i_pa._arguments.size(); i++)
		{
			const _ProgramArgument &a = i_pa._arguments[i];
			std::cout << "'" << a.key << "' => '" << a.value << "' (processed=" << a.argumentParsedAndAccessed << ")" << std::endl;
		}
		return io_os;
	}
};

}

#endif
