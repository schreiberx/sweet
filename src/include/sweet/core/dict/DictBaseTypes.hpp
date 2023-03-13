/*
 * Dict.hpp
 *
 *  Created on: Feb 18, 2023
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_DICT_BASE_HPP__
#define SRC_INCLUDE_SWEET_DICT_BASE_HPP__

#include <string>
#include <vector>
#include <complex>
#include <fstream>
#include <sweet/core/SWEETError.hpp>
#include "DictArrayND.hpp"

namespace sweet
{

/**
 * A class with information only available to Dict implementations
 */
class DictBaseTypes
{
public:
	typedef long long int64;
	typedef double float64;
	typedef std::complex<double> complex128;
};
}	// namespace

#endif
