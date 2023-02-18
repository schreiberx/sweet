/*
 * Dict.hpp
 *
 *  Created on: Feb 13, 2023
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_DICT_HPP__
#define SRC_INCLUDE_SWEET_DICT_HPP__

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/SWEETError.hpp>
#include "DictArrayND.hpp"
#include "DictBaseTypes.hpp"
#include "DictElements.hpp"
#include "DictFileReadWrite.hpp"

namespace sweet
{

/**
 * A simple dictionary which allows storing scalar and array data
 *
 * It also supports reading/writing to/from files in C++ and Python
 */
class Dict: public DictBaseTypes
{
private:
	std::vector<_DictElement> _dict;

public:
	int size()	const
	{
		return _dict.size();
	}

public:
	const _DictElement& operator[](int i_index)	const
	{
		return _dict[i_index];
	}
	_DictElement& operator[](int i_index)
	{
		return _dict[i_index];
	}

private:
	bool _debug;

	const char* SWEET_MAGIC_FILE_CODE = "SWEETFileDict";


private:
	void _basicTest()
	{
		if (sizeof(int64) != 8)
			SWEETError("Something is weird with int64!");

		if (sizeof(float64) != 8)
			SWEETError("Something is weird with float64!");
	}

public:
	Dict(bool i_debug = false)	:
		_debug(i_debug)
	{
		_basicTest();
	}

	Dict(const std::string &i_filename, bool i_debug = false)
	{
		_debug = i_debug;

		_basicTest();

		fileLoad(i_filename);
	}

public:
	template <typename T>
	void get(const std::string &i_key, T &o_value)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
			{
				_dict[i].get(o_value);
				return;
			}
		}

		std::ostringstream ss;
		ss << "Key '" << i_key << "' not found!";
		SWEETError(ss.str());
	}



	bool keyExists(const std::string &i_key)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return true;
		}

		return false;
	}

	int keyIndex(const std::string &i_key)	const
	{
		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			if (_dict[i].getKey() == i_key)
				return i;
		}

		return -1;
	}


	template <typename T>
	void set(const std::string &i_key, const T &i_value, bool i_ignore_existing = true)
	{
		int idx = keyIndex(i_key);

		if (idx < 0)
		{
			_dict.push_back(_DictElement(i_key, i_value));
			return;
		}

		if (!i_ignore_existing)
			SWEETError("Key already exists");

		_dict[idx].set(i_key, i_value);
	}


public:
	friend
	std::ostream& operator<<(
			std::ostream& io_os,
			const Dict &i_dict
	)
	{
		for (std::size_t i = 0; i < i_dict._dict.size(); i++)
		{
			const _DictElement &e = i_dict._dict[i];

			std::cout << " + " << e << std::endl;
		}

		return io_os;
	}


	/*
	 * Load from a Dict
	 *
	 * See corresponding Python file for description of file format
	 */
public:
	bool fileLoad(const std::string &i_filename)
	{
		DictFileRead f(i_filename);

		std::string magic_start = f.loadStr0();

		if (magic_start != SWEET_MAGIC_FILE_CODE)
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_start << "' at beginning!";
			SWEETError(ss.str());
		}

		std::size_t num_entries = f.loadData<int64>();


		if (_debug)
			std::cout << "Dict: Found " << num_entries << " dictionary entries" << std::endl;

		_dict.resize(num_entries);
		for (std::size_t i = 0; i < num_entries; i++)
		{
			_DictElement &e = _dict[i];

			e.fileLoadKeyTypeValue(f);

			if (_debug)
				std::cout << " + Loaded '" << e.getKey() << "' => '" << e.getValueAsString() << "'" << std::endl;
		}

		std::string magic_end = f.loadStr0();

		if (magic_start != SWEET_MAGIC_FILE_CODE)
		{
			std::ostringstream ss;
			ss << "Invalid Magic code '" << magic_end << "' at end!";
			SWEETError(ss.str());
		}

		return false;
	}


	/*
	 * Write to a Dict
	 */
public:
	bool fileSave(const std::string &i_filename, int i_verbosity = 0)
	{
		DictFileWrite f(i_filename);

		f.writeStr0(SWEET_MAGIC_FILE_CODE);

		f.writeData<int64>(_dict.size());

		for (std::size_t i = 0; i < _dict.size(); i++)
		{
			_DictElement &e = _dict[i];

			e.fileWriteKeyTypeValue(f);

			if (_debug)
				std::cout << " + Loaded '" << e.getKey() << "' => '" << e.getValueAsString() << "'" << std::endl;
		}

		f.writeStr0(SWEET_MAGIC_FILE_CODE);

		return false;
	}

};

}	// namespace

#endif
