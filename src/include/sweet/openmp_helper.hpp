/*
 * openmp_helper.hpp
 *
 *  Created on: 20 Jul 2015
 *      Author: Martin Schreiber <schreiberx@gmail.com>
 */
#ifndef SRC_EXAMPLES_OPENMP_HELPER_HPP_
#define SRC_EXAMPLES_OPENMP_HELPER_HPP_


#if !SWEET_SIMD_ENABLE

	#define OPENMP_SIMD

#else

	#ifndef __GNUC__
		#define OPENMP_SIMD	simd
	#else
		#if __GNUC__ > 4 || ((__GNUC__ == 4) && (__GNUC_MINOR__ > 8))
			#define OPENMP_SIMD	simd
		#else
			#warning "SIMD is disabled for this compiler version"
			#define OPENMP_SIMD
		#endif
	#endif

#endif

#endif /* SRC_EXAMPLES_OPENMP_HELPER_HPP_ */
