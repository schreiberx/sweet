/*
 * openmp_helper.hpp
 *
 *  Created on: 20 Jul 2015
 *      Author: Martin Schreiber <SchreiberX@gmail.com>
 */
#ifndef SRC_EXAMPLES_OPENMP_HELPER_HPP_
#define SRC_EXAMPLES_OPENMP_HELPER_HPP_


/**
 * This is a class to overcome SIMD instruction issues with older GNU compilers
 */

#define OMP_SCHEDULE	schedule(static)

#if !SWEET_THREADING

	#define SWEET_OMP_PAR_FOR_SIMD

#else

	#if !SWEET_SIMD_ENABLE

		#define OPENMP_PAR_SIMD	OMP_SCHEDULE
		#define OPENMP_SIMD
		#define SWEET_OMP_PAR_FOR_SIMD #pragma omp parallel

	#else

		#ifdef __INTEL_COMPILER
			#define OPENMP_PAR_SIMD     simd OMP_SCHEDULE
			#define OPENMP_SIMD simd
		#else
			#define OPENMP_PAR_SIMD	simd OMP_SCHEDULE
			#define OPENMP_SIMD simd
		#endif


		/*
		 * This should be used for all parallel for loops with stream-like access
		 * e.g.
		 * #pragma omp OMP_PAR_FOR_SIMD
		 */

		#define OMP_PAR_FOR_SIMD	parallel for OPENMP_PAR_SIMD
		#define SWEET_OMP_PAR_FOR_SIMD #pragma omp parallel for OPENMP_PAR_SIMD
	#endif

#endif



#endif /* SRC_EXAMPLES_OPENMP_HELPER_HPP_ */
