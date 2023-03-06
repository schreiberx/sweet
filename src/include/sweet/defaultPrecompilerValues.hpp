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
 *
 * We set them to default values here to support the IDE.
 */

#ifndef SWEET_GUI
	#error "SWEET_GUI not set"
	#define SWEET_GUI 1
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


#endif
