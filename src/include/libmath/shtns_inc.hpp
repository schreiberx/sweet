/*
 * shtns_inc.hpp
 *
 *  Created on: 25 Aug 2016
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_
#define SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_


//#define SHTNS_COMPLEX_SPH_OLD_INTERFACE	1
#define SHTNS_COMPLEX_SPH_OLD_INTERFACE	0

/**
 * This is a wrapper to avoid including shtns.h0 twice.
 * #ifndef stuff is missing in the shtns.h0 header file
 */

// undef for C++11 support
// TODO: For which compiler versions is this required?
#undef _COMPLEX_H
#include <shtns.h>


#endif /* SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_ */
