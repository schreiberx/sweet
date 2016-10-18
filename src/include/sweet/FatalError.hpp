/*
 * FatalError.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: martin
 */

#ifndef SRC_INCLUDE_SWEET_FATALERROR_HPP_
#define SRC_INCLUDE_SWEET_FATALERROR_HPP_


class FatalError
{
public:
	FatalError(const std::string i_error)
	{
		std::cerr << i_error << std::endl;
		assert(false);
		exit(-1);
	}
};


#endif /* SRC_INCLUDE_SWEET_FATALERROR_HPP_ */
