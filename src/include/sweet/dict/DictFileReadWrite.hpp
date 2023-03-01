/*
 * Dict.hpp
 *
 *  Created on: Feb 18, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_DICT_FILE_HPP__
#define SRC_INCLUDE_SWEET_DICT_FILE_HPP__

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/SWEETError.hpp>
#include "DictArrayND.hpp"
#include "DictBaseTypes.hpp"

namespace sweet
{

/**
 * A class implementing file-based convenient functions
 */
class DictFileRead
{
	std::ifstream is;

public:
	DictFileRead(const std::string &i_filename)
	{
		is = std::ifstream(i_filename, std::ios::in | std::ios::binary);

		if (!is.is_open())
			SWEETError(std::string("Unable to open file ")+i_filename);
	}


public:
	/**
	 * Read 0 terminated string from file
	 */
	std::string loadStr0()
	{
		std::vector<char> buffer;
		buffer.resize(1024);

		bool found = false;
		for (std::size_t i = 0; i < buffer.size()-1; i++)
		{
			is.read((char*)&(buffer[i]), sizeof(char));
			if (buffer[i] == 0)
			{
				found = true;
				break;
			}
		}

		if (!found)
			SWEETError("0 termination character not found, stopping here");

		// convert to std::string
		return (std::string)(char*)buffer.data();
	}

	template <typename T>
	T loadData()
	{
		T retval;
		is.read((char*)&retval, sizeof(retval));
		return retval;
	}

	template <int D, typename T>
	void loadArray(
			sweet::DictArrayND<D,T> &array
	)
	{
		for (std::size_t i = 0; i < array.size(); i++)
			array.data()[i] = loadData<T>();
	}
};


/**
 * A class implementing file-based convenient functions
 */
class DictFileWrite
{
	std::ofstream os;

public:
	DictFileWrite(const std::string &i_filename)
	{
		os = std::ofstream(i_filename, std::ios::out | std::ios::binary);

		if (!os.is_open())
			SWEETError(std::string("Unable to open file ")+i_filename);
	}


public:
	/**
	 * Write 0 terminated string from file
	 */
	void writeStr0(const std::string& i_value)
	{
		os.write(i_value.c_str(), i_value.length()+1);
	}

	template <typename T>
	void writeData(const T& i_value)
	{
		os.write((const char*)&i_value, sizeof(i_value));
	}

	template <int D, typename T>
	void writeDictArrayRawData(
			const sweet::DictArrayND<D,T> &i_array
	)
	{
		for (std::size_t i = 0; i < i_array.size(); i++)
			writeData<T>(i_array.data()[i]);
	}

	template <int D, typename T>
	void writeDictArray(
			const sweet::DictArrayND<D,T> &i_array
	)
	{
		// Write out shape
		for (int d = 0; d < D; d++)
			writeData<DictBaseTypes::int64>(i_array.shape()[d]);

		writeDictArrayRawData(i_array);
	}
};

}	// namespace

#endif
