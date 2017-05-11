/*
 * FatalError.hpp
 *
 *  Created on: 18 Oct 2016
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_SWEET_FATALERROR_HPP_
#define SRC_INCLUDE_SWEET_FATALERROR_HPP_


#include <cassert>
#include <string>
#include <iostream>

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



class AssertFatalError
{
public:
#ifdef NDEBUG
	AssertFatalError(bool i_assertion, const std::string i_error)	{}
#elif !SWEET_DEBUG
	AssertFatalError(bool i_assertion, const std::string i_error)	{}
#else
	AssertFatalError(bool i_assertion, const std::string i_error)
	{
		if (i_assertion)
			return;

		std::cerr << std::flush << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << "ASSERT ERROR: " << i_error << std::endl;
		std::cerr << "********************************************" << std::endl;
		std::cerr << std::endl;
		assert(false);
		exit(-1);
	}
#endif

};


#endif /* SRC_INCLUDE_SWEET_FATALERROR_HPP_ */
