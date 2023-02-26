/*
 * StringSplit.hpp
 *
 *  Created on: 12 Feb 2017
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_
#define SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_

#include <vector>
#include <string>
#include <string>

#include "SWEETError.hpp"

class StringSplit
{
public:
	static
	std::vector<std::string> split(
			const std::string &i_string,
			const std::string &i_delimiter
	)
	{
		std::vector<std::string> ret;

		// create a copy of the string
		std::string s = i_string;

		// http://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
		std::size_t pos = 0;
		while ((pos = s.find(i_delimiter)) != std::string::npos) {
			std::string found_string = s.substr(0, pos);
		    ret.push_back(found_string);
		    s.erase(0, pos + i_delimiter.length());
		}

		ret.push_back(s);

		return ret;
	}

public:
	/**
	 * @brief fill a vector of doubles with a string that contains substrings separated by a delimiter
	 * 
	 * If the string contains a different number of elements than o_doubles, function throws an error.
	 * Exception: If the string contains one element, o_doubles gets filled up with the value.
	 * 
	 * @param i_str : the input string (example: "2,3,0.1")
	 * @param o_doubles : the output vector to be filled (example: here it should have 3 elements)
	 * @param i_delimiter : the delimiter that should be used (example: ",")
	 */
	static
	void split_n_doubles(
			const std::string & i_str,
			std::vector<double> & o_doubles,
			const std::string & i_delimiter
	)
	{
		std::vector<std::string> substrings = split(i_str, i_delimiter);

		if ((substrings.size() != 1) && (substrings.size() != o_doubles.size()))
		{
			SWEETError("This number of values is unexpected! String '" + i_str + "' invalid.");
		}

		if (substrings.size() == 1)
		{
			// fill o_doubles with the given value
			for (auto & value : o_doubles)
			{
				value = atof(i_str.c_str());
			}
			return;
		}
		// both vectors have the same length
        for (size_t iter = 0; iter < o_doubles.size(); iter++)
        {
            o_doubles.at(iter) = atof(substrings.at(iter).c_str());
        }
	}

	/**
	 * @brief fill a vector of integers with a string that contains substrings separated by a delimiter
	 * 
	 * If the string contains a different number of elements than o_ints, function throws an error.
	 * Exception: If the string contains one element, o_ints gets filled up with the value.
	 * 
	 * @param i_str : the input string (example: "2,3,1")
	 * @param o_ints : the output vector to be filled (example: here it should have 3 elements)
	 * @param i_delimiter : the delimiter that should be used (example: ",")
	 */
	static
	void split_n_ints(
			const std::string & i_str,
			std::vector<int> & o_ints,
			const std::string & i_delimiter
	)
	{
		std::vector<std::string> substrings = split(i_str, i_delimiter);

		if ((substrings.size() != 1) && (substrings.size() != o_ints.size()))
		{
			SWEETError("This number of values is unexpected! String '" + i_str + "' invalid.");
		}

		if (substrings.size() == 1)
		{
			// fill o_ints with the given value
			for (auto & value : o_ints)
			{
				value = atoi(i_str.c_str());
			}
			return;
		}
		// both vectors have the same length
        for (size_t iter = 0; iter < o_ints.size(); iter++)
        {
            o_ints.at(iter) = atoi(substrings.at(iter).c_str());
        }
	}


	static
	int split3double(
			const std::string &i_str,
			double *o_int0,
			double *o_int1,
			double *o_int2
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atof(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			return 2;
		}
		else if (c == 3)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			*o_int2 = atof(res[2].c_str());
			return 3;
		}

		SWEETError("More than 3 values given");
		return -1;
	}


	static
	int split2int(
			const std::string &i_str,
			int *o_int0,
			int *o_int1
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atoi(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atoi(res[0].c_str());
			*o_int1 = atoi(res[1].c_str());
			return 2;
		}

		SWEETError("More than 2 values given");
		return -1;
	}


	static
	int split2double(
			const std::string &i_str,
			double *o_int0,
			double *o_int1
	)
	{
		std::vector<std::string> res = StringSplit::split(i_str, ",");

		int c = res.size();

		if (c == 0)
			SWEETError("Invalid format for modes");

		if (c == 1)
		{
			*o_int0 = atof(res[0].c_str());
			return 1;
		}
		else if (c == 2)
		{
			*o_int0 = atof(res[0].c_str());
			*o_int1 = atof(res[1].c_str());
			return 2;
		}

		SWEETError("More than 2 values given");
		return -1;
	}


};



#endif
