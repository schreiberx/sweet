#ifndef CGL_SHADER_TEXTURIZE_RAINBOW_MAP_HPP
#define CGL_SHADER_TEXTURIZE_RAINBOW_MAP_HPP

#include <libgl/core/GlError.hpp>
#include <libgl/core/GlTexture.hpp>


class GlShaderTexturizeRainbowmap	: public GlProgram
{
public:
	GlUniform pvm_matrix_uniform;

	GlShaderTexturizeRainbowmap()
	{
		if (!initVertFragShadersFromDirectory("shader_texturize_rainbowmap"))
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

	~GlShaderTexturizeRainbowmap()
	{
	}
};


#endif
