/*
 * shtns_inc.hpp
 *
 *  Created on: 25 Aug 2016
 *      Author: Martin Schreiber <SchreiberX@gmail.com> Schreiber <M.Schreiber@exeter.ac.uk>
 */

#ifndef SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_
#define SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_


/**
 * This is a wrapper to avoid including shtns.h twice.
 * #ifndef stuff is missing in the shtns.h header file
 */

// undef for C++11 support
// TODO: For which compiler versions is this required?
#undef _COMPLEX_H
#include <shtns.h>


#endif /* SRC_INCLUDE_LIBMATH_SHTNS_INC_HPP_ */
