#ifndef CGL_SHADER_COLOR_HPP
#define CGL_SHADER_COLOR_HPP

#include "libgl/core/CGlError.hpp"

/**
 * shader for a single color without anything else
 */
#include "libgl/core/CGlProgram.hpp"


class GlShaderColor	:
	public CGlProgram
{
	CGlShader &getVertShaderSingleton()
	{
		static CGlShader vertShader;
		return vertShader;
	}

	CGlShader &getFragShaderSingleton()
	{
		static CGlShader fragShader;
		return fragShader;
	}

	bool &getLoadedSingleton()
	{
		static bool loaded = false;
		return loaded;
	}

	int &getUsageCounterSingleton()
	{
		static int usage_counter = 0;
		return usage_counter;
	}

public:
	CGlUniform frag_color;	///< uniform to vec4 specifying fragment color
	CGlUniform pvm_matrix;	///< uniform to proj-view-model matrix


	GlShaderColor()
	{
		std::string infoLog;

		if (!getLoadedSingleton())
		{
			getVertShaderSingleton().init(GL_VERTEX_SHADER);

			if (!getVertShaderSingleton().loadSource(SHADER_GLSL_DEFAULT_DIR"shader_color/vertex_shader.glsl"))
				return;

			if (!getVertShaderSingleton().compile())
				return;

			getFragShaderSingleton().init(GL_FRAGMENT_SHADER);

			if (!getFragShaderSingleton().loadSource(SHADER_GLSL_DEFAULT_DIR"shader_color/fragment_shader.glsl"))
				return;

			if (!getFragShaderSingleton().compile())
				return;

			getLoadedSingleton() = true;
		}

		this->attachShader(getVertShaderSingleton());
		this->attachShader(getFragShaderSingleton());

		// link programs
		if (!link())
		{
			std::string infoLog;
			getInfoLog(infoLog);
			std::cerr << "info Log: during linking: " << infoLog << std::endl;
			return;
		}

		bindAttribLocation(0, "vertex_position");
		setupUniform(frag_color, "frag_color");
		setupUniform(pvm_matrix, "pvm_matrix");

		getUsageCounterSingleton()++;
	}

	~GlShaderColor()
	{
		getUsageCounterSingleton()--;

		if (getUsageCounterSingleton() == 0)
		{
			getLoadedSingleton() = false;
			getVertShaderSingleton().freeIfValid();
			getFragShaderSingleton().freeIfValid();
		}
	}
};
#endif
