/*
 * StringSplit.hpp
 *
 *  Created on: 12 Feb 2017
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_
#define SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_

#include <vector>
#include <string>
#include <string>

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
};



#endif /* SRC_INCLUDE_SWEET_STRINGSPLIT_HPP_ */
