/*
 * defaultPrecompilerValues.hpp
 *
 *  Created on: Feb 19, 2023
 *      Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_INCLUDE_SWEET_DEFAULTPRECOMPILERVALUES_HPP_
#define SRC_INCLUDE_SWEET_DEFAULTPRECOMPILERVALUES_HPP_

/*
 * These precompiler directives should be automatically set by the compile command.
 * If not, they will trigger an error.
 *
 * We also use this to set them to default values here to support the IDE.
 */

#ifndef SWEET_GUI
	#error "SWEET_GUI not set"
	#define SWEET_GUI 1
#endif

#ifndef SWEET_EIGEN
	#error "SWEET_EIGEN not set"
	#define SWEET_EIGEN 1
#endif

#ifndef SWEET_PARAREAL
	#error "SWEET_PARAREAL not set"
	#define SWEET_PARAREAL 1
#endif

#ifndef SWEET_LIBPFASST
	#error "SWEET_LIBPFASST not set"
	#define SWEET_LIBPFASST 1
#endif

#ifndef SWEET_XBRAID
	#error "SWEET_XBRAID not set"
	#define SWEET_XBRAID 1
#endif

#ifndef SWEET_USE_SPHERE_SPECTRAL_SPACE
	#error "SWEET_USE_SPHERE_SPECTRAL_SPACE not set"
	#define SWEET_USE_SPHERE_SPECTRAL_SPACE 1
#endif

#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "SWEET_USE_PLANE_SPECTRAL_SPACE not set"
	#define SWEET_USE_PLANE_SPECTRAL_SPACE 1
#endif



/*
 * This precompiler directive should !!! NEVER !!! be used
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_SPACE
	#define SWEET_USE_PLANE_SPECTRAL_SPACE	1
#endif

#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#if SWEET_USE_PLANE_SPECTRAL_DEALIASING
		#error "Dealiasing only available with plane-spectral-space enabled"
	#endif
#endif

/*
 * Activating this option via the precompiler also activates the FFTW library
 */
#ifndef SWEET_USE_LIBFFT
	#define SWEET_USE_LIBFFT	1
#endif

/*
 * Activating the dealiasing creates plans and additional information for dealiasing strategies
 */
#ifndef SWEET_USE_PLANE_SPECTRAL_DEALIASING
	#define SWEET_USE_PLANE_SPECTRAL_DEALIASING 1
#endif



#endif
