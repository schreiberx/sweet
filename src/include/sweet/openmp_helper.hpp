/*
 * openmp_helper.hpp
 *
 *  Created on: 20 Jul 2015
 *      Author: Martin Schreiber <M.Schreiber@exeter.ac.uk>
 */
#ifndef SRC_EXAMPLES_OPENMP_HELPER_HPP_
#define SRC_EXAMPLES_OPENMP_HELPER_HPP_


/**
 * This is a class to overcome SIMD instruction issues with older GNU compilers
 */

#define OMP_SCHEDULE	schedule(static)

#if !SWEET_SIMD_ENABLE

	#define OPENMP_PAR_SIMD	OMP_SCHEDULE
	#define OPENMP_SIMD

#else

	#ifdef __INTEL_COMPILER
		#define OPENMP_PAR_SIMD     simd OMP_SCHEDULE
		#define OPENMP_SIMD simd
	#else
		#define OPENMP_PAR_SIMD	simd OMP_SCHEDULE
		#define OPENMP_SIMD simd
	#endif
#endif


#endif /* SRC_EXAMPLES_OPENMP_HELPER_HPP_ */
