/**
 * math include wrapper to avoid problems for missing definitions of M_PIl on apple machines
 */

#include <cmath>

// apple hack
#ifndef M_PIl
#	define M_PIl M_PI
#endif
