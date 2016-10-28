/*
 * FatalError.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_FATALERROR_HPP_
#define SRC_INCLUDE_SWEET_FATALERROR_HPP_


#include <cassert>

class FatalError
{
public:
	FatalError(const std::string i_error)
	{
		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << "ERROR: " << i_error << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		assert(false);
		exit(-1);
	}
};


#endif /* SRC_INCLUDE_SWEET_FATALERROR_HPP_ */
