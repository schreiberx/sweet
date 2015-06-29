#ifndef CGL_SHADER_TEXTURIZE_HPP
#define CGL_SHADER_TEXTURIZE_HPP

#include "libgl/core/CGlTexture.hpp"
#include "libgl/core/CGlError.hpp"


class GlShaderTexturize	: public CGlProgram
{
public:
	CGlUniform pvm_matrix_uniform;

	GlShaderTexturize()
	{
		if (!initVertFragShadersFromDirectory("shader_texturize"))
			return;

		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: linking: " << infoLog << std::endl;
			return;
		}

		setupUniform(pvm_matrix_uniform, "pvm_matrix");
	}

	~GlShaderTexturize()
	{
	}
};


#endif
