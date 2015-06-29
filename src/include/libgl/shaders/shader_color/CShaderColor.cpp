#include "CShaderColor.hpp"


/**
 * if this file is included by multiple cpp files which are compiled, the following lines have to be moved to a separate cpp file
 * to be linked only once! otherwise a compiler error would be generated
 */
#if 0
// no shaders are initialized so far
CGlShader GlShaderColor::vertShader;
CGlShader GlShaderColor::fragShader;
bool GlShaderColor::loaded = false;
int GlShaderColor::usage_counter = 0;
#endif
