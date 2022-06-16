/*
 * StringSplit.hpp
 *
 *  Created on: 12 Feb 2017
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
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
	 * Exception: If the string contains one element, o_double gets filled up with the value.
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
};



#endif /* SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_ */
